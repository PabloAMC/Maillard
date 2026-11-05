# Bork, Haase, Rohn & Kanzler 2022 (Food Chem. 380:131852) — EXTRACTION (paper key `bork2022`; the PDF on disk is `data/articles/bork2022b.pdf`, see sec. 0)
### *"Formation of melanoidins – Aldol reactions of heterocyclic and short-chain Maillard intermediates"*
### Norfuraneol (NF) 20 mM and methylglyoxal (MGO) or diacetyl (DA) ~20 mM, each alone or 1:1, in **unbuffered water adjusted to pH 5.0 ± 0.1** (NaOH / HCl), sealed ampoules, **130 °C, 0-300 min**; NF by HPLC-DAD (external calibration), MGO as its OPD quinoxaline; colour index (A420 per mmol); SEC; HRMS of the oligomers. **Fig. 1D is a six-point time course of NF heated ALONE at pH 5 / 130 °C: first-order, k ≈ 4.8-5.1 × 10⁻³ min⁻¹ (t½ ≈ 140 min) over 0-300 min** — the first fed-NF decay at acid pH in the corpus, from a laboratory not yet in the repo.

**Read 2026-11-05 (first reader) from `data/articles/bork2022b.pdf`** (10 pp., Elsevier publisher PDF, clean text layer).
Prose from `pdftotext -raw`; **Fig. 1C / 1D digitised from a 400 dpi crop of p. 4** (scratchpad `img/bork22b_fig1d-04.png`,
`img/bork22b_fig1c-04.png`; marker centres located by colour, script `scripts/digit.py`). Read error **±1.5 % of the
100 mol % scale** (marker size and axis calibration). Fits: scratchpad `scripts/bork_fit.py` (maillard env). Supplementary
material (Fig. S-1 colour of sugar / amino acid pots; **Fig. S-2 NF / DA in the NF/DA pot**; Figs. S-3 to S-8; Tables
S-1 to S-3 HRMS assignments) is **NOT on disk**.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Leon V. Bork, Paul T. Haase, Sascha Rohn, **Clemens Kanzler** (corresponding, clemens.kanzler@tu-berlin.de) |
| Affiliation | Technische Universität Berlin, Institute of Food Technology and Food Chemistry, Dept. of Food Chemistry and Analysis, Gustav-Meyer-Allee 25, 13355 Berlin |
| Venue | *Food Chemistry* **380** (2022) 131852 |
| DOI | **10.1016/j.foodchem.2021.131852** (printed p. 1) — matches the brief |
| Dates | received 31 Aug 2021; revised 2 Dec 2021; accepted 9 Dec 2021; online 14 Dec 2021 |

[!] **File naming is crossed.** The brief listed `bork2022.pdf` as this paper; the file of that name is the OTHER Bork 2022
paper (Food Chem. 395:133592, NF + glyoxal / glycolaldehyde / acetaldehyde; dossier `bork2022b_extraction.md`). This paper
is in `bork2022b.pdf`. Paper keys follow the brief and the house order (the earlier paper, vol. 380, online Dec 2021, is
`bork2022`; vol. 395, online Jun 2022, is `bork2022b`). `data/keys/pdf_dossier_aliases.yml` maps the two file stems
crosswise. Renaming the two PDFs in the main checkout would let the alias lines go.

**Laboratory.** Kanzler's group, TU Berlin (Food Chemistry and Analysis; Rohn chair). **Not in
`data/keys/laboratories.yml`**; proposed id **`tu_berlin_kanzler`**, lab string "TU Berlin (Kanzler / Rohn)". Same group
as `bork2022b` and `bork2023`; **count the three papers as ONE laboratory**. Not the BAM group of `mohsin2018` (Kanzler a
co-author there; `docs/reference/contrast_benchmark/lab_ids.yml` keeps "TU Berlin / BAM (Kanzler)" for that one) — the
conservative leave-one-lab-out choice would merge them; the evidence (Mohsin 2018 is a BAM materials-analysis paper) says
apart. Owner call.

## 1. What this paper answers for the model

| model need | what Bork 2022 gives | tag |
|---|---|---|
| NF decay vs time at pH 5, 130 °C, a lab not yet in the repo | **Fig. 1D, NF alone, 6 points 0-300 min**: 100 / 85.9 / 75.1 / 59.2 / 42.2 / 20.0 mol % | [F] |
| first-order k and window | **k = 4.8 × 10⁻³ min⁻¹ (weighted, ±0.2 stat.) to 5.1 × 10⁻³ (unweighted), 0-300 min**; pointwise 4.4-5.4 × 10⁻³, no trend over 0-180 min; the 300 min point sits ~2σ below the line (weak sign of late acceleration) | [Z] |
| buffered or drifting pH | **Unbuffered**; pH set to 5.0 ± 0.1 at room temperature before heating; **no pH measured after heating** | [M] / [!] |
| O2 status | Sealed ampoules; **no degassing, no chelator, headspace and ampoule volume not stated** (air presumed) | [!] |
| MGO alone, pH 5 / 130 °C | Fig. 1D: 100 / 33.8 / 22.0 / 13.7 / 10.4 / 7.3 mol % — **not first order** (pointwise k 0.036 → 0.0087 min⁻¹) | [F] / [Z] |
| NF + MGO cross reaction | Fig. 1C, 1:1: NF 19.2 / 12.0 / 7.5 / 4.7 / 1.9 %, MGO 4.0 / 1.9 / 2.3 / 1.6 / 1.4 % at 30-300 min | [F] |
| NF + DA | text only: DA ~40 % left at 300 min, NF "fully converted" (Fig. S-2, not on disk) | [M] / [Q] |
| NF fragments / products | HRMS: aldol addition / condensation oligomers of NF, MGO or DA, with redox steps (±H2); NF's oxidised triketone form invoked (Kanzler 2016 [C]) | [Q] / [P] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots (sec. 2.2) | NF **57 mg** (0.50 mmol) and MGO 40 wt% **90 µL** (~0.58 mmol [Z]: density ~1.17 g/mL assumed) or DA **43 µL** (0.49 mmol [Z]) in **25 mL water** → **NF 20 mM**, MGO ~23 mM, DA ~20 mM [Z]; each alone or as the binary mixture |
| pH | "after adjusting the pH to 5.0 ± 0.1" [V]; NaOH / HCl (the companion paper names 1 M of each); **label = initial, at room temperature; no buffer; not re-measured** |
| Vessel / heating | "sealed ampoules", 130 °C (heating medium not stated); times 0, 30, 60, 120, 180, 300 min; quench in a freezer at −20 °C |
| Headspace / O2 | not stated; no N2 / Ar; no DTPA / EDTA in the pot (EDTA only in the HPLC eluent) |
| Low-water arm (sec. 2.7) | same masses, 60 wt% water (i.e. 40 wt% solids), 130 °C, 0-120 min; colour only |
| Replicates | triplicate pots, mean ± SD |

## 3. Assays [M]

- **NF:** HPLC-DAD 285 nm, Nucleosil C18, ion-pair eluent (5 mM phosphate pH 6.0, 12.5 mM EDTA, 2.5 mM TBAHS)/MeOH;
  **external calibration** (method of Kanzler 2017 [C]). Results plotted as mol % of the initial amount.
- **MGO:** OPD 50 mM, 1:1 with sample, 24 h dark RT; quinoxaline by HPLC-DAD 318 nm, external calibration (Bruhns 2018 [C]).
- [!] Both retention times are printed as "tR = 41.9 min" — one is a typo; immaterial to the numbers.
- **No LOD / LOQ** for any analyte. Error bars only visible at 120 min (NF ±3 %) and 30 min (MGO ±6 %); elsewhere smaller
  than the markers (≤ ~2 %).
- Colour index = A420 × dilution / Σ mmol of reactants.

## 4. Data

### 4.1 Fig. 1D — individual incubations, water pH 5.0 (initial), 130 °C [F, ±1.5 mol %]

| t (min) | 0 | 30 | 60 | 120 | 180 | 300 |
|---|---|---|---|---|---|---|
| **NF alone (mol %)** | 100 | 85.9 | 75.1 | 59.2 (±3 SD) | 42.2 | 20.0 |
| MGO alone (mol %) | 100 | 33.8 (±6 SD) | 22.0 | 13.7 | 10.4 | 7.3 |

### 4.2 Fig. 1C — NF + MGO 1:1, same conditions [F, ±1.5 mol %]

| t (min) | 0 | 30 | 60 | 120 | 180 | 300 |
|---|---|---|---|---|---|---|
| NF (mol %) | 100 | 19.2 | 12.0 | 7.5 | 4.7 | 1.9 |
| MGO (mol %) | 100 | 4.0 | 1.9 | 2.3 | 1.6 | 1.4 |

Text: after MGO is "nearly consumed at around 30 min", NF had fallen by more than 80 % [M]. The authors read the similar
conversion of both as an equimolar reaction [Q].

### 4.3 Arithmetic [Z]

1. **NF first-order decay, alone, pH 5 (initial), 130 °C.** ln-linear fit through 100 % at t = 0, weighted by read error:
   **k = 4.80 × 10⁻³ ± 0.18 × 10⁻³ min⁻¹** (t½ 145 min); free intercept 4.83 × 10⁻³ (C0 = 100.4 %); unweighted 5.1 × 10⁻³.
   Pointwise ln(C0/C)/t: 5.07 / 4.77 / 4.37 / 4.79 / 5.36 × 10⁻³ at 30 / 60 / 120 / 180 / 300 min. **First order holds
   to read error over 0-180 min**; the 300 min residual is −0.17 in ln (≈ 2σ of read error at 20 %), a weak hint of
   acceleration late (cf. Nakamura 2020's accelerating NF decay at pH 8, `nakamura2020_extraction.md` sec. 4a).
   **Quoted value: k(130 °C, pH 5 initial, unbuffered, air presumed) = 4.9 (± 0.4) × 10⁻³ min⁻¹, window 0-300 min**
   (the ± covers weighting choice and read error, not lab-to-lab spread; see sec. 6 for the same lab's repeat).
2. **No come-up correction**: ampoule heating time not stated; with a 25 mL fill the lag could be several minutes. If the
   first ~5 min are cold, k rises by ≤ ~2 % (5/300) for the long window and ~17 % for the 30-min point. Immaterial to the
   order of magnitude.
3. **MGO alone is not first order**: pointwise k 0.036 → 0.025 → 0.017 → 0.013 → 0.0087 min⁻¹. A single k is not a
   description; ~66 % lost in 30 min, then a slow tail (hydrate / oligomer equilibria or a second-order sink [P]).
4. **NF + MGO cross reaction, crude bound.** Treating the extra NF loss in the 1:1 pot over 0-30 min as bimolecular with
   equal concentrations: 1/C − 1/C0 = (1/3.84 − 1/20) mM⁻¹ = 0.21 mM⁻¹ in 30 min → k2 ≳ **7 × 10⁻³ mM⁻¹ min⁻¹** at
   130 °C, pH 5 (lower bound: MGO is already 96 % gone by 30 min, so the true early rate is higher). Not a fit row.

## 5. What it can and cannot test in the model

**Can:**
- **The engine's NF sink at acid pH.** The shipped `r_nf_decay` is pH-free, log10 k(145 °C) = −2.182 at the carbonyl-sink
  barrier 143.5 kJ/mol (`results/validation/core_prediction_uncertainty.md` row `b8.k_nf_decay`;
  `src/kinetic_core/parameters_sulfur.py` family B) → **1.42 × 10⁻³ min⁻¹ at 130 °C**. Bork measures **4.9 × 10⁻³: ×3.4
  the engine** [Z]. B61's hydroxide term (`NF_HYDROXIDE_DECAY`, OFF; 0.032 min⁻¹ at in-situ pH 8.08 / 100 °C, a decade per
  pH unit, declared 103.7 kJ/mol) gives **3.2 × 10⁻⁴ min⁻¹ at pH 5 / 130 °C — 16× below Bork**; engine pH-free + B61
  together: 1.7 × 10⁻³, still ×2.9 below. **At pH 5 the loss is not hydroxide-dominated**; a pH-independent (or weakly
  pH-dependent) channel carries it.
- **The pH slope between Mikami (pH 8, 100 °C) and Bork (pH 5, 130 °C).** Moving Bork to 100 °C needs a barrier, which
  nobody has measured: with 81.5 / 103.7 / 143.5 kJ/mol, k(100 °C, pH 5) = 7.2 / 4.2 / 1.6 × 10⁻⁴ min⁻¹, so Mikami's
  0.032 is 44 / 75 / 196× faster over 3 pH units = **0.55-0.76 decade per pH unit** [Z] — sub-decade, i.e. the
  specific-base term of B61 plus a parallel pH-flat term, which is what the engine's two NF sinks already are in form;
  Bork sets the size of the flat one. Caveats: two labs, two temperatures, buffered (Mikami) vs unbuffered (Bork),
  1.7 vs 20 mM.
- **Lau 2003's bound** (k_d ≲ 6 × 10⁻³ min⁻¹ at 2σ, whey gel, 116 °C, pH n.r.). Bork moved to 116 °C: 2.1 / 1.7 / 1.1 ×
  10⁻³ min⁻¹ (81.5 / 103.7 / 143.5 kJ/mol) — **consistent with Lau** for any of these barriers, and also with Pandit 2006's
  "≤ ~6 × 10⁻³". Neither gel bound tests Bork.
- **MGO's own loss at pH 5 / 130 °C** (Fig. 1D) and **an NF + MGO aldol sink** (Fig. 1C; the engine has no NF + MGO step:
  `docs/reference/pathway_atlas/F31_melanoidins.yml` scope, not checked line by line).

**Cannot:**
- **No temperature series** — no Ea for NF decay. One pH (5.0, start only). One NF concentration (20 mM): self-reaction
  (NF + NF aldol, NF + its own MGO / acetaldehyde / formaldehyde fragments — the companion paper's HRMS shows all of them)
  cannot be separated from a unimolecular loss, although first-order kinetics over a five-fold fall argues against a
  dominant second-order term.
- **pH drift unknown.** In unbuffered water NF breakdown likely makes acids (formic / acetic [P]); the pot is probably
  below pH 5 after the first hour, which would only push the hydroxide term further down.
- **O2 unknown.** Air-saturated water holds ~0.26 mM O2 at 25 °C [Z, handbook value] = 6.5 µmol in 25 mL, against
  ~400 µmol NF lost by 300 min: dissolved O2 alone can carry ≤ ~2 % of the loss at 1:1. A headspace of air of the same
  volume would add ~210 µmol O2 (≤ ~50 %). So an oxidative channel (Nakamura 2020 Fig. 7 shows DTPA and N2 suppress NF
  loss at pH 8) **cannot be excluded or confirmed**; the constant k over 300 min argues against an O2-limited channel
  unless the headspace is large. The authors themselves invoke oxidation of NF to a triketone [Q].
- Products of NF alone are not quantified (no MGO, acetol, formaldehyde, acids measured in the NF-alone pot).

## 6. Conflicts and agreements with the repo's records

1. **[!] CONTRADICTS the engine's NF sink at pH 5** (×3.4 over the shipped pH-free `r_nf_decay`; ×16 over B61's hydroxide
   term alone). Same direction as Mikami vs engine at pH 8 (×700, `mikami2017_extraction.md`) but much smaller. The
   `engine_vs_literature.yml` row F13 cites only Mikami / Reading for NF decay; Bork is a third, acid-side point.
   Not registered as a CF here (`docs/reference/` is owned by another agent today); proposed CF text: "NF decay at pH 5 /
   130 °C, Bork 2022 Fig. 1D, 4.9e-3 min⁻¹ vs engine 1.4e-3 (pH-free) + 3e-4 (B61); one lab, unbuffered, air".
2. **SUPPORTS B61's structure, not its sufficiency**: a specific-base term alone cannot carry NF loss at pH 5.
3. **Same-lab repeat, `bork2023_extraction.md`**: NF alone, same nominal pot (20 mM, pH 5.0, 130 °C, ampoules), **48 ± 4 %
   left at 120 min vs 59.2 ± 3 % here** → k(0-120) 6.1 × 10⁻³ vs 4.4 × 10⁻³ min⁻¹: a **within-lab, between-paper spread of
   ~0.2 in ln (×1.4 in k)**, larger than either paper's replicate SD. Bork 2023 used a different HPLC method and
   quantified relative to t = 0. MGO alone at 120 min: 16 % (2023) vs 13.7 % (here) — agree.
4. **Same-lab, possibly the same data, `bork2022b_extraction.md`**: its text prints NF 19 % at 30 min in NF + MGO, equal
   to Fig. 1C here (19.2 %) [!] duplicate_suspected — enter once.
5. **Hofmann 1998 / Hofmann-Schieberle NF levels** are net (formation − decay) in amine pots; this paper's fed-NF decay is
   the first acid-side sink measurement to set against them; no direct conflict.

## 7. Proposed ledger entries (trunk.yml; NOT written — first reader only, second reader pending)

```yaml
  - id: bork2022_fig1d_norfuraneol_alone_decay_ph5_130c
    dossier: bork2022_extraction.md
    paper: "Bork, Haase, Rohn & Kanzler 2022 Food Chem 380:131852"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "norfuraneol 20 mM alone (57 mg in 25 mL)"
    medium: "water, no buffer, pH adjusted to 5.0 +- 0.1 with NaOH/HCl before heating; sealed ampoules; headspace not stated"
    temperature_c: [130]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [0, 30, 60, 120, 180, 300]
    species: [norfuraneol]
    kind: time_series
    location: "Fig. 1D p. 4 (digitised in dossier sec. 4.1, +-1.5 mol %)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: true
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "heated separately or as mixtures in sealed ampoules at 130 C after adjusting the pH to 5.0"
    flags: [ph_start_only, unbuffered, o2_not_stated, no_lod, same_lab_repeat_differs_bork2023]
    note: "First-order k 4.8-5.1e-3 1/min over 0-300 min (dossier 4.3). x3.4 engine pH-free r_nf_decay at 130 C; x16 B61 hydroxide term at pH 5."
  - id: bork2022_fig1d_methylglyoxal_alone_ph5_130c
    dossier: bork2022_extraction.md
    paper: "Bork, Haase, Rohn & Kanzler 2022 Food Chem 380:131852"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "methylglyoxal ~23 mM alone (90 uL of 40 wt% in 25 mL)"
    medium: "as above"
    temperature_c: [130]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [0, 30, 60, 120, 180, 300]
    species: [methylglyoxal]
    kind: time_series
    location: "Fig. 1D p. 4 (dossier sec. 4.1)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: true
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "Methylglyoxal was quantified as its quinoxaline derivative ... by external calibration"
    flags: [ph_start_only, unbuffered, o2_not_stated, not_first_order]
    note: "66 % lost in 30 min, then a slow tail; MGO concentration from the 40 wt% label (density assumed)."
  - id: bork2022_fig1c_norfuraneol_methylglyoxal_1to1_ph5_130c
    dossier: bork2022_extraction.md
    paper: "Bork, Haase, Rohn & Kanzler 2022 Food Chem 380:131852"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "norfuraneol 20 mM + methylglyoxal ~23 mM"
    medium: "as above"
    temperature_c: [130]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [0, 30, 60, 120, 180, 300]
    species: [norfuraneol, methylglyoxal]
    kind: time_series
    location: "Fig. 1C p. 4 (dossier sec. 4.2)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "After methylglyoxal was nearly consumed at around 30 min"
    flags: [ph_start_only, unbuffered, o2_not_stated, duplicate_suspected_bork2022b_text]
    note: "Engine has no NF + MGO step; k2 >~ 7e-3 /mM/min (dossier 4.3 item 4). bork2022b prints the same 19 % at 30 min."
```

## 8. Supplementary material

Not on disk. Wanted: **Fig. S-2** (NF and DA in the NF/DA pot vs time; the NF-alone curve may be replotted there) and
Tables S-1 to S-3 (HRMS assignments). Nothing else model-relevant is missing.
