# Bork, Haase, Rohn & Kanzler 2022b (Food Chem. 395:133592) — EXTRACTION (paper key `bork2022b`; the PDF on disk is `data/articles/bork2022.pdf`, see sec. 0)
### *"Structural characterization of polar melanoidins deriving from Maillard reaction intermediates – A model approach"*
### Norfuraneol (NF) 20 mM alone or 1:1 with glyoxal (GO), glycolaldehyde (GA), acetaldehyde (AA) or furfural (FF) in **unbuffered water adjusted to pH 5** (1 M NaOH / HCl), sealed ampoules, **130 °C, 30-300 min**; low-water arm (40 wt% solids); colour index, SEC, HRMS / Kendrick analysis of the oligomers; NF's proposed degradation scheme (Fig. 3). **Every NF time course is in Fig. S-2 (Supplementary; read 2026-11-05, sec. 7)**; the main text prints only the NF left at 30 min in each binary pot. **The SI's NF-alone curve is `bork2022`'s Fig. 1D replotted (with larger error bars): this paper adds no independent NF-alone decay.**

**Read 2026-11-05 (first reader) from `data/articles/bork2022.pdf`** (pp. 1-10, Elsevier publisher PDF, clean text layer).
Prose from `pdftotext -raw`. No figure in the main text carries an NF concentration, so nothing was digitised.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Leon V. Bork, Paul T. Haase, Sascha Rohn, **Clemens Kanzler** (corresponding) — TU Berlin, Food Chemistry and Analysis |
| Venue | *Food Chemistry* **395** (2022) 133592 |
| DOI | **10.1016/j.foodchem.2022.133592** (printed p. 1) — matches the brief's second Bork 2022 |
| Dates | received 30 Mar 2022; revised 20 Jun 2022; accepted 26 Jun 2022; online 28 Jun 2022 |

[!] **The file `bork2022.pdf` is THIS paper, not Food Chem. 380:131852** (that one is `bork2022b.pdf`, dossier
`bork2022_extraction.md`). Keys follow the brief and publication order; the crossing is recorded in
`data/keys/pdf_dossier_aliases.yml`.

**Laboratory.** Same group as `bork2022` and `bork2023` (proposed `tu_berlin_kanzler`); one laboratory.

## 1. What this paper answers for the model

| model need | what Bork 2022b gives | tag |
|---|---|---|
| NF decay alone vs time, pH 5, 130 °C | **Not in the main text; SI Fig. S-2 replots `bork2022` Fig. 1D** (85.5 / 74.9 / 59.0 / 42.1 / 20.2 % at 30-300 min; SD ±6 at 30, ±4.7 at 120 min) — sec. 7.2 | [F] / [!] duplicate |
| NF left at 30 min with a partner (1:1, 20 mM each) | **GO 0.7 %**, **MGO 19 %**, **GA 60 %**, **FF 80 %**, **AA 83 %** of the initial NF | [M] |
| reactivity ranking in water (colour index + NF conversion) | GO > MGO > AA > GA > FF | [Q] |
| NF degradation products | proposed: hydration / ring opening to a 2,4-diketone → retro-aldol to **MGO + GA**; oxidative hydrolysis to a **triketone** → retro-aldol to **formaldehyde + acetol**; HRMS of dry-heated NF (125 °C, 10 min) shows NF oligomers with MGO, acetaldehyde and formaldehyde units and oxidised NF | [P] / [Q] |
| pH, O2 | as `bork2022`: unbuffered, pH set to 5 before heating; no degassing, no chelator stated | [!] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots (sec. 2.2) | NF 57 mg (0.50 mmol) with GO 40 wt% 72.5 µL, GA 30 mg (0.50 mmol), AA 22 mg (0.50 mmol) or FF 48.0 mg (0.50 mmol) in **25 mL water** → **20 mM each** [Z]; alone or binary |
| pH | "at pH 5" [V]; 1 M NaOH / 1 M HCl; no buffer; **initial, room temperature; not re-measured** |
| Heating | sealed ampoules, 130 °C, 30 / 60 / 120 / 180 / 300 min; freezer quench −20 °C; centrifuged, diluted |
| Low-water arm (sec. 2.6) | same masses at 60 wt% water, 130 °C, 10 / 30 / 60 / 120 min; colour of water and ethyl acetate extracts |
| Replicates | triplicate, mean ± SD |
| NF assay (sec. 2.4) | HPLC-DAD 285 nm, ion-pair (phosphate pH 6 / EDTA / TBAHS) on Nucleosil C18; **external calibration and relative to t = 0** |
| LOD / LOQ | **not given** |

## 3. Data [M]

- Text, sec. 3.1: NF left after **30 min at 130 °C, pH 5**: "0.7 % with regard to the initial concentration" with GO;
  "19 % (methylglyoxal), 60 % (glycolaldehyde), 80 % (furfural), and 83 % (acetaldehyde)" [V, abridged].
- Fig. 1 (colour index after 120 min, aqueous pH 5 and 40 wt% water): bars, not digitised (colour, not a model observable
  here).
- SEC (Fig. 2, 300 min): NF/GO colorants span > 800 kDa to < 0.34 kDa; NF with AA / GA / MGO / DA limited to ~1-10 kDa
  and ~0.34 kDa fractions [Q].

**Arithmetic [Z].** NF alone at 30 min in the same nominal pot is **85.9 %** (`bork2022` Fig. 1D). So the partners add an
apparent first-order NF loss over 0-30 min of ln(0.859/x)/30: **GO 0.16, MGO 0.050, GA 0.012, FF 2.4 × 10⁻³, AA 1.2 ×
10⁻³ min⁻¹** (FF and AA at the read-error limit of the NF-alone value). These are pot-specific (20 mM partner) and not
rate constants of a step.

## 4. What it can and cannot test

**Can:** rank C2 / C3 carbonyls as NF sinks at pH 5 / 130 °C (GO ≫ MGO ≫ GA > FF ≈ AA); supply the authors' product
scheme for NF degradation (MGO + GA by retro-aldol; formaldehyde + acetol via an oxidised triketone) — mechanism [P],
consistent with Mikami 2017's MGO / diacetyl from NF at pH 8 only in the MGO part (diacetyl is not in this scheme).

**Cannot:** give independent NF-alone kinetics (the SI's is `bork2022`'s, sec. 7.2); give fragment yields (none measured); separate oxidative from
hydrolytic NF loss (no O2 control).

## 5. Conflicts and agreements with the repo's records

1. [!] **duplicate, confirmed by the SI**: "19 % (methylglyoxal)" at 30 min equals `bork2022` Fig. 1C (19.2 %), and SI
   Fig. S-2D equals Fig. 1C at every time (sec. 7.2). The same incubation reported twice; never enter both.
2. **Engine:** the NF sinks are `r_nf_decay` (→ inert FRAG_C), B61's `r_nf_decay_oh` (+ optional MGO share) and the H2S
   routes; there is no NF + GO / GA / MGO aldol sink. The fast NF + GO loss (t½ ≈ 4 min at 20 mM each [Z]) is a candidate
   missing sink wherever GO and NF co-exist (pentose pots), not tested here.
3. The proposed retro-aldol NF → MGO + GA agrees in direction with B61's `mgo_share` structure (`parameters_sulfur.py`
   sub-model 3); no yield is given to set it.

## 6. Proposed ledger entry (trunk.yml; NOT written — first reader only, second reader pending)

```yaml
  - id: bork2022b_text_norfuraneol_left_30min_with_carbonyls_ph5_130c
    dossier: bork2022b_extraction.md
    paper: "Bork, Haase, Rohn & Kanzler 2022 Food Chem 395:133592"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "norfuraneol 20 mM + glyoxal, methylglyoxal, glycolaldehyde, furfural or acetaldehyde 20 mM (1:1)"
    medium: "water, no buffer, pH 5 set with 1 M NaOH/HCl before heating; sealed ampoules"
    temperature_c: [130]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [30]
    species: [norfuraneol]
    kind: endpoint
    location: "text sec. 3.1 p. 3; full time courses in SI Fig. S-2 (dossier sec. 7.2)"
    status: text_only
    quant_class: calibrated
    simulable: false
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "only traces of norfuraneol (0.7 % with regard to the initial concentration)"
    flags: [ph_start_only, unbuffered, o2_not_stated, single_time_point, duplicate_confirmed_bork2022_fig1c_via_si]
    note: "NF left at 30 min: GO 0.7, MGO 19, GA 60, FF 80, AA 83 %. NF alone in the same pot: bork2022 Fig. 1D (85.9 % at 30 min)."
  - id: bork2022b_si_figs2_norfuraneol_with_carbonyls_time_course_ph5_130c
    dossier: bork2022b_extraction.md
    paper: "Bork, Haase, Rohn & Kanzler 2022 Food Chem 395:133592 (Supplementary Fig. S-2)"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "norfuraneol 20 mM + glyoxal, glycolaldehyde, acetaldehyde or furfural 20 mM (1:1)"
    medium: "water, no buffer, pH 5 set with 1 M NaOH/HCl before heating; sealed ampoules; headspace not stated"
    temperature_c: [130]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [0, 30, 60, 120, 180, 300]
    species: [norfuraneol]
    kind: time_series
    location: "SI Fig. S-2 A-C, E (digitised in dossier sec. 7.2, +-1 mol %)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "individual and combined heat treatment of norfuraneol and different carbonyl compounds"
    flags: [ph_start_only, unbuffered, o2_not_stated, partner_not_measured, nf_alone_curve_is_bork2022_fig1d]
    note: "Panel D (MGO) is bork2022 Fig. 1C and is NOT entered here. Engine has no NF + carbonyl aldol sink. NF + GA first order (k ~1.9e-2), NF + AA accelerating, NF + FF ~6.7e-3 1/min."
```

## 7. Supplementary material (read 2026-11-05, first reader)

**Source.** `data/articles/bork2022_supplementary.docx` (Word file, Figs. S-1 to S-16, Tables S-1 to S-5). **The file
named `bork2022_` belongs to THIS paper (key `bork2022b`)**, crossed like the PDFs (sec. 0). Identified by content: its
header carries this paper's exact title and authors, and its figures are this paper's pots (NF with GO / GA / AA / MGO /
FF, dry-heated NF, NF / GO and NF / GA at 40 °C). Read with `textutil` and a zip / XML reader (tables); Fig. S-2 is an
embedded TIFF (1323 × 1940 px, five panels), flattened on white and digitised by a square-marker finder (two dark
vertical edges + top and bottom edges on each time column; filled vs open by interior grey), axes from the tick marks
(2.86-2.87 px per mol %, 1.67 px per min); error bars from their caps. **Read error ±1 mol %.** Scripts: scratchpad
`bork_tools/squares.py`, `errbars.py`, `si_fit.py` (maillard env, `-I`).

### 7.1 Items asked for

| item | found? | value |
|---|---|---|
| ampoule / vial volume, fill, headspace, degassing, atmosphere, chelator / metals | **NOT found** | no methods text in the SI |
| pH after heating | **NOT found** | captions: "pH 5" (initial) only |
| **NF alone vs time** | **FOUND — a replot, not new data** | dashed open squares in every panel of Fig. S-2: 85.5 / 74.9 / 59.0 / 42.1 / 20.2 mol % at 30 / 60 / 120 / 180 / 300 min, 130 °C, pH 5 (initial) = `bork2022` Fig. 1D within 0.6 mol %. **New: the replicate SDs, ±6 at 30 min and ±4.7 at 120 min** (others hidden by the marker, ≤ 2.5) |
| NF temperature or pH series | **NOT found** | 130 °C and pH 5 only (plus 40 °C HRMS pots and a dry 125 °C pot, no kinetics) |

[!] Captions give the pots as **"0.5 wt%"**; the methods' masses (57 mg NF in 25 mL; partner 0.50 mmol) are 0.23 wt%
NF and ~0.35-0.45 wt% total solids [Z]. Taken as a rounded label; the 20 mM of sec. 2 stands.

### 7.2 Fig. S-2 — NF vs time, 20 mM NF ± 20 mM partner, water pH 5 (initial), 130 °C [F ±1 mol %]

Caption (paraphrased): relative NF concentration during individual and combined heating of NF with different carbonyls
(0.5 wt%, pH 5) at 130 °C; solid = combined, dashed = NF alone. SD from error-bar caps where visible.

| t (min) | 0 | 30 | 60 | 120 | 180 | 300 |
|---|---|---|---|---|---|---|
| NF alone (dashed, every panel) | 100 | 85.5 (SD ≈ 6) | 74.9 | 59.0 (SD ≈ 4.7) | 42.1 | 20.2 |
| (A) NF + glyoxal | 100 | **0.7** | 0.2 | 0.2 | 0.2 | 0.2 |
| (B) NF + glycolaldehyde | 100 | 59.4 (SD ≈ 7) | 31.1 | 10.0 | 3.7 | 0.3 |
| (C) NF + acetaldehyde | 100 | 83.1 (SD ≈ 6) | 67.0 | 38.4 (SD ≈ 4.5) | 20.8 | 5.6 |
| (D) NF + methylglyoxal | 100 | 19.2 (SD ≈ 3) | 12.2 | 7.5 | 4.7 | 1.6 |
| (E) NF + furfural | 100 | 80.8 | 68.4 | 45.2 | 35.6 (SD ≈ 5) | 12.0 |

The 30 min values match the text of sec. 3 (0.7 / 60 / 83 / 19 / 80 %). Panel D equals `bork2022` Fig. 1C point for
point (19.2 / 12.0 / 7.5 / 4.7 / 1.9): **the duplicate of sec. 5 item 1 is confirmed** — one incubation, two papers.

**Arithmetic [Z]** (pointwise ln(100/C)/t, × 10⁻³ min⁻¹, at 30 / 60 / 120 / 180 / 300 min):
1. **NF alone**: 5.2 / 4.8 / 4.4 / 4.8 / 5.3 → **k = 4.91 × 10⁻³ ± 0.23 × 10⁻³ min⁻¹** (SD-weighted, through 100 %;
   χ²/dof 0.40; window 0-300 min). Same data as `bork2022` sec. 4.3 (4.8-5.1 × 10⁻³); **no conflict, not independent**.
   The 300 min point sits ~1σ below the line with these SDs (it was ~2σ against read error alone).
2. **NF + GA: first order in NF throughout**, 17 / 20 / 19 / 18 / 19 → k ≈ 1.9 × 10⁻² min⁻¹, i.e. +1.4 × 10⁻² over NF
   alone, constant for 300 min although GA is equimolar. Either GA is regenerated or the sink is pseudo-first-order
   in NF through a species that does not run down [P]; GA itself was not measured.
3. **NF + AA accelerates**: 6.2 / 6.7 / 8.0 / 8.7 / 9.6 — a product-driven (autocatalytic-looking) sink [P].
4. **NF + FF: first order**, 7.1 / 6.3 / 6.6 / 5.7 / 7.1 → k ≈ 6.7 × 10⁻³ (+1.8 × 10⁻³ over alone).
5. **NF + MGO / GO**: front-loaded; most loss within 30 min, then the slow tail (MGO, GO spent).
6. The sec. 3 increments over 0-30 min, recomputed against the SI's NF-alone 85.5 %: GO 0.16, MGO 0.050, GA 0.012,
   FF 1.9 × 10⁻³, AA 0.9 × 10⁻³ min⁻¹ — unchanged in substance.

### 7.3 Every figure and table

| item | caption (paraphrased) | content / values | tag |
|---|---|---|---|
| Fig. S-1 | colour index of the partners heated ALONE, 130 °C, 120 min: (A) water (0.5 wt%, pH 5), (B) 40 wt% water, aqueous vs ethyl acetate extracts | read ±0.01 / ±0.005 mmol⁻¹: (A) GO ≈ 0 (SD 0.11), GA 0, **FF 1.56 (SD 0.20)**, AA 0.23 (SD 0.76); (B) aqueous / EA: GO 0.021 / 0.001, GA 0.315 (SD 0.05) / 0.019, FF 0.279 (SD 0.10) / 0.234 (SD 0.14), AA 0.023 / 0.002. Colour only | [F] |
| **Fig. S-2** | NF vs time alone and with GO / GA / AA / MGO / FF | sec. 7.2 | [F] |
| Fig. S-3 | ESI(+)-HRMS of NF heated dry, 125 °C, 10 min: scan and MS² of m/z 355, 363, 379, 509 | spectra | [Q] |
| Figs. S-4 to S-7 | proposed fragmentation pathways of m/z 355, 363, 379, 509 | schemes | [P] |
| Fig. S-8 | HRMS of NF / GA in water, 40 °C, 100 min: scan and Kendrick plots (H2O, NF, GA bases) | plots | [Q] |
| Figs. S-9 / S-10 | Kendrick plots of NF / GO and NF / GA normalised to NF fragments (MGO, oxidised NF, acetaldehyde, formaldehyde) | plots; fragment units present in the oligomers | [Q] |
| Fig. S-11 | Kendrick plots of the NF-alone HRMS spectrum (H2O, NF, MGO, oxidised NF, acetaldehyde, formaldehyde) | plots | [Q] |
| Figs. S-12 to S-15 | MS² / MS³ of NF / GO ions (m/z 309, 367, 425, 499, 557, 595, 653) and NF / GA m/z 197, 40 °C, 100 min | spectra | [Q] |
| Fig. S-16 | UHPLC-DAD (285 nm) and single-ion chromatograms (m/z 197, 257) of NF / GA, 40 °C, 100 min, with MS² | chromatograms | [Q] |
| **Table S-1** | HRMS of dry-heated NF, 125 °C, 10 min, 100 most intense signals | 100 ions, each assigned as n NF + MGO / acetaldehyde / formaldehyde units ± H2O ± H2. **99 of 100 ions (98 % of summed intensity) carry a net loss of H2** (oxidised); base peak 3 NF + MGO + 3 H2O − 3 H2 (m/z 485, Na+); next 2 NF + MGO + 2 H2O − 2 H2 (355), 3 NF + H2O − 2 H2 (379), 4 NF + 2 H2O − 3 H2 (509) | [M] / [Z] count |
| Table S-2 | assignment of the MS² spectra of Fig. S-3 | fragment ions of m/z 355, 363, 379, 509 | [M] |
| Table S-3 | HRMS of NF / GO, 40 °C, 100 min | 100 ions; base peak 2 NF + GO (m/z 309); 40 ions (18 % of intensity) oxidised | [M] / [Z] count |
| Table S-4 | HRMS of NF / GA, 40 °C, 100 min | 100 ions; base peak NF + GA (m/z 197); 9 ions (1 % of intensity) oxidised | [M] / [Z] count |
| Table S-5 | MS² / MS³ assignments for Figs. S-14 / S-15 | fragment ions | [M] |

**What the SI adds for the model.** (i) The replicate SDs of the NF-alone series (±6 at 30, ±4.7 at 120 min).
(ii) Five NF-sink time courses with partner carbonyls at pH 5 / 130 °C — shapes differ (first order with GA and FF,
accelerating with AA, front-loaded with GO / MGO), so the 30-min increments of sec. 3 are not rate constants.
(iii) **Indirect on O2:** in NF heated dry (125 °C, 10 min, atmosphere not stated), nearly every oligomer is
dehydrogenated (−H2), and the dominant ones contain an MGO unit — NF's own self-reaction runs through oxidised species
and its retro-aldol MGO [Q]. Whether the oxidant is O2 or internal redox (disproportionation) is not tested; and it is a
dry pot, not the pH 5 solution. **Nothing on O2 status, final pH or temperature dependence of NF-alone loss.**
