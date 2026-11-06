# Kanzler, Falkenhagen, Rüger, Mier & Rohn 2026 (Food Chem. 529:151175) — EXTRACTION (file stem `kanzler2026`)
### *"Amino acids promote melanoidin formation through aldol reactions in the MAILLARD chemistry of methylglyoxal"*
### Methylglyoxal (MGO) **100 mM**, alone or 1:1 with one of ten amino acids (Ala, Val, Ile, Phe, Pro, Ser, Glu, Gln, **Lys**, Arg; **no glycine, no cysteine**), in **unbuffered water adjusted to pH 5.0 ± 0.1**, 2.5 mL in sealed glass ampoules, **100 °C, 0-120 min**; plus a 60-140 °C series at 30 min (Ala, Phe, Pro), an Ala series 10-100 mM, and a water-content series (SI). MGO by OPD / HPLC-DAD with **external calibration**; amino acids by HPTLC relative to t = 0. **MGO alone: 65 % left at 30 min, 31 % at 120 min** — consistent with second order in MGO (k2 ≈ 1.85 × 10⁻⁴ mM⁻¹ min⁻¹), first-order k falling 0.014 → 0.008 min⁻¹. Amino acids speed MGO loss 2-60× at 30 min while most are barely consumed. **No norfuraneol (NF).**

**Read 2026-11-05 (first reader) from `data/articles/kanzler2026.pdf`** (13 pp., Elsevier open-access PDF, clean text
layer). Prose from `pdftotext -raw`; **Figs. 1D, 1E, 2B, 2D, 2F digitised from 300 dpi crops** of pp. 4 and 6 (scratchpad
`sweep_img/k26f1d-04.png`, `k26f1e-04.png`, `k26f2b-06.png`, `k26f2d-06.png`, `k26f2f-06.png`; marker centres and bar tops
read against the axis ticks, 4.05-4.24 px per %): **read error ±1.5 % of the 100 % scale**. Other panels from 130 dpi
renders (±3 %). Arithmetic: scratchpad `tools/sweep_fits.py` (maillard env, `-I`). Supplementary material
(`kanzler2026_supplementary.docx`) read into **sec. 8**.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[D]** read from the vector drawing embedded in the SI (exact coordinates; equivalent to [F] at ~0.5 % precision) ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **Clemens Kanzler** (corresponding, clemens.kanzler@tu-berlin.de), Johann Falkenhagen, Tatjana Rüger, Leon Valentin Mier, Sascha Rohn |
| Affiliation | TU Berlin, Institute of Food Technology and Food Chemistry, Dept. of Food Chemistry and Analysis, Kaiserin-Augusta-Allee 14, 10553 Berlin |
| Venue | *Food Chemistry* **529** (2026) 151175 |
| DOI | **10.1016/j.foodchem.2026.151175** (printed p. 1) — matches the brief |
| Dates | received 20 Jan 2026; revised 10 Sep 2026; accepted 12 Sep 2026; online 14 Sep 2026 |

**Laboratory.** Kanzler group, TU Berlin (Food Chemistry and Analysis; new address, same department) — **the same
laboratory as `bork2022`, `bork2022b`, `bork2023`, `bork2024`, `bork2024b`** (proposed `tu_berlin_kanzler`); count
them as one laboratory for leave-one-lab-out.

## 1. What this paper answers for the model

| model need | what Kanzler 2026 gives | tag |
|---|---|---|
| MGO loss alone, pH 5, 100 °C | **100 / 65 (SD ≈ 10) / 31 (SD ≈ 6.5) % at 0 / 30 / 120 min**, 100 mM; text: "around 35% of MGO were consumed" at 30 min, "around 31%" remained at 120 min | [F] / [M] |
| order of that loss | two points: k1 0.0144 (0-30) → 0.0082 min⁻¹ (30-120); k2 1.79 / 1.87 × 10⁻⁴ mM⁻¹ min⁻¹ — **second order fits, first order does not; order not identified within the SD** | [Z] |
| MGO loss with each amino acid | Fig. 1D (30 min) and 1E (120 min) for all ten; full time courses (0 / 30 / 60 / 120) for Ala, Phe, Pro | [F] |
| amine dependence | MGO + Ala 10-100 mM: excess k over MGO alone ∝ [Ala]^0.75-0.84 | [Z] |
| temperature | 30 min at 60-140 °C, MGO + Ala / Phe / Pro: apparent Ea(80 → 100 °C) 121 / 95 / 97 kJ/mol | [Z] |
| amino acids consumed | 120 min: Ala ~0 %, Lys ~2 %, Glu ~10 %, Phe ~37 %, Arg ~52 %, Gln ~85 % (to pyroglutamate) | [F] / [M] |
| glycine, cysteine, NF | **none** | — |
| products | ESI(−)-HRMS: aldol / Michael oligomers of MGO and its fragments formaldehyde, acetaldehyde, acetol, glycolaldehyde, pyruvic acid; MG-H1 / MG-H3 with Arg; pyroglutamate with Gln | [M] / [P] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots (sec. 2.2) | MGO 40 wt% (Merck) and amino acid, **each 100 mmol/L, in water**; MGO alone "under otherwise identical conditions" |
| pH | "adjusted to (5.0 ± 0.1)" with 1 M HCl / NaOH [V, paraphrased]; **label = initial, at room temperature; no buffer; not re-measured after heating** |
| Vessel / heating | **2.5 mL** sealed in a glass ampoule (ampoule volume not stated); heating block at (100 ± 1) °C; 0, 30, 60, 120 min; cooled in a −20 °C freezer |
| Temperature series | Ala, Phe, Pro pots, 30 min at (60 / 80 / 120 / 140 ± 1) °C (plus 100 °C) |
| Concentration series | Ala 10, 25, 50 mmol/L with MGO 100 mM (plus 100 mM), 30 / 60 / 120 min, 100 °C |
| Water-content series | MGO / Ala, 40-98.4 wt% water, 100 °C, 30 min (sec. 3.3 and SI Fig. S3; preparation not described in the methods) [!] |
| Headspace / O2 / metals | **not stated**: no degassing, no inert gas, no chelator; air presumed |
| Replicates | triplicate pots, mean ± SD; one- or two-way ANOVA + Tukey |

## 3. Assays [M]

- **MGO:** 25 µL sample + 475 µL water + 500 µL OPD (50 mM, water / MeOH 1:1), overnight dark RT; 2-methylquinoxaline by
  HPLC-DAD 318 nm (ProntoSil Phenyl), tR 42.1 min. **External calibration with 2-methylquinoxaline**, standards
  equivalent to **1-150 mmol/L in the pot** (diluted 1:40 like the samples). No LOD / LOQ printed.
  - [!] The **lowest calibrant is 1 % of C0**: every MGO value below ~1 % (Pro at 60-120 min, Glu at 120 min, the 120 /
    140 °C points) is extrapolated below the calibrated range; treat < 1 % as "≤ 1 %".
  - [!] Calibrating with the quinoxaline assumes complete derivatisation of the (hydrated / oligomeric) MGO.
- **Amino acids:** HPTLC on cellulose, ninhydrin, 570 nm (400 nm for Pro), **relative to the unheated sample on the same
  plate**. No LOD.
- Colour A420 × dilution (aqueous; ethyl acetate extract at 120 min); CIELab; TEAC; SEC (PL aquagel-OH, 420 nm);
  ESI(−)-Orbitrap HRMS.
- **Replicate noise** (SD from error bars): MGO alone ±10 (30 min) / ±6.5 (120 min) % of C0 — large; MGO with amino
  acids ±0.5-4.5 %; Ala series ±0.5-1.5 %.

## 4. Data

### 4.1 Fig. 1D / 1E — MGO remaining (% of C0) after 30 / 120 min, 100 mM each, pH 5 (initial), 100 °C [F ±1.5 %; SD]

| partner | − (alone) | Ala | Val | Ile | Phe | Pro | Ser | Glu | Gln | **Lys** | Arg |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 min | **65 (10)** | 32 (1.5) | 31.5 (2) | 27.5 (1) | 22 (2) | 1.5 (0.5) | 18.5 (1) | 4 (1) | 19 (0.5) | **11 (4.5)** | 16 (2.5) |
| 120 min | **31 (6.5)** | 4.3 (0.5) | 3.7 (0.5) | 2.8 (0.5) | 1.6 (0.5) | ≤ 1 | 2.2 (0.5) | ≤ 1 | 2.6 (0.5) | **5.6 (2.5)** | 8.4 (1.8) |

Text: MGO conversion "at least 68%" after 30 min and "over 90%" at 120 min in all amino-acid pots [M] (agrees: the
highest 30 min value is 32 %).

### 4.2 Fig. 1F — amino acid remaining after 120 min (% of t = 0) [F ±3 %; SD]

Ala 100 (1), Val 84 (7), Ile 89.5 (8), Phe 63 (7), Pro 79 (2), Ser 87 (30), Glu 90 (3.5), Gln 15 (2), **Lys 98 (5)**,
Arg 47.5 (0.5). Text: Lys "around 2%" converted; Arg "more than 50%"; Glu 10 %; Gln to ~15 %; Phe "converted to 36%"
[M] — [!] read as 36 % *consumed* (the figure shows 63 % left).

### 4.3 Fig. 2B — MGO remaining vs time, 100 mM each, 100 °C [F ±1.5 %]

| t (min) | 0 | 30 | 60 | 120 |
|---|---|---|---|---|
| MGO + Ala | 100 | 31.8 | 14.0 | 4.3 |
| MGO + Phe | 100 | 21.6 | 6.2 | 1.7 |
| MGO + Pro | 100 | 1.9 | ≤ 1 | ≤ 1 |

### 4.4 Fig. 2D — MGO remaining after 30 min vs temperature, 100 mM each [F ±1.5 %; SD where visible]

| T (°C) | 60 | 80 | 100 | 120 | 140 |
|---|---|---|---|---|---|
| MGO + Ala | 99.5 | 86.5 (1) | 26.5 | 2.5 | ≤ 1 |
| MGO + Phe | 93 (4) | 72 (5) | 15.5 | ≤ 1 | ≤ 1 |
| MGO + Pro | 96 | 53.5 | 2.5 | ≤ 1 | ≤ 1 |

Text: Pro at 80 °C "reduced to 53%" [M]. [!] **The 100 °C / 30 min points differ from Fig. 1D / 2B**: Ala 26.5 vs 31.8,
Phe 15.5 vs 21.6 — a separate incubation set; between-set spread ~5-6 % of C0 (~0.2 in ln at this level), larger than
the within-set SDs.

### 4.5 Fig. 2F — MGO remaining vs Ala concentration, MGO 100 mM, 100 °C [F ±1.5 %]

| Ala (mM) | 10 | 25 | 50 | 100 |
|---|---|---|---|---|
| 30 min | 57.2 (1.3) | 49.0 (1) | 40.4 | 31.8 |
| 60 min | 38.9 | 30.4 | 22.1 | 14.3 |
| 120 min | 23.1 | 16.3 | 9.6 | 4.2 |

[!] Text says "(F) MGO conversion after 120 min"; the figure and caption show 30, 60 and 120 min. The 100 mM column equals
Fig. 2B within read error (same pots or a repeat).

### 4.6 Colour, TEAC [F ±3 % of panel scale; colour only]

Fig. 1A, A420 × F after 120 min: alone 1.0; Ala 13.3; Val 11.2; Ile 10.1; Phe 6.3; **Pro 32**; Ser 14.1; Glu 14.4; Gln 14.5;
**Lys 22.5**; Arg 3.6. Fig. 1C, TEAC (mmol TE/L): alone 7; amino-acid pots 21-29; Arg 35. No engine observable.

### 4.7 Arithmetic [Z] (`sweep_fits.py`)

1. **MGO alone, 100 mM, pH 5 (initial), 100 °C.** First order through the two points: k1(0-30) = 0.0144, k1(30-120) =
   0.0082, k1(0-120) = 0.0098 min⁻¹ — decelerating. Second order: k2 = 1.79 × 10⁻⁴ (0-30) and 1.87 × 10⁻⁴ mM⁻¹ min⁻¹
   (0-120) — **constant**. Two-point apparent order n = 1.9 nominal; **at ±1 SD of the 30 / 120 min points (±10 / ±6.5)
   n spans ~0.6 to unbounded**, so the order is suggested, not identified. Quoted value: **k2(MGO + MGO) ≈ 1.85 × 10⁻⁴
   mM⁻¹ min⁻¹ at 100 °C, pH 5 initial, 100 mM, window 0-120 min** (equivalently pseudo-first-order 0.010-0.014 min⁻¹ at
   100 mM).
2. **Same lab, other temperature:** `bork2022` Fig. 1D MGO alone (23 mM, pH 5 initial, 130 °C): pointwise k2 2.84 → 1.84 ×
   10⁻³ mM⁻¹ min⁻¹ over 30-300 min (also decelerating as second order). Against 1.85 × 10⁻⁴ at 100 °C: **apparent Ea
   ≈ 95-113 kJ/mol** for the MGO self-loss — two papers, two concentrations, two fill volumes (25 vs 2.5 mL); an
   indication, not a fit.
3. **Amino acids at 100 mM, k1(0-30):** Ala 0.038, Val 0.038, Ile 0.043, Phe 0.051, Pro 0.14, Ser 0.056, Glu 0.11, Gln
   0.055, Lys 0.074, Arg 0.061 min⁻¹ vs 0.0144 alone: **×2.6 (Ala) to ×9.7 (Pro)**, Lys ×5. With Ala the interval k1 falls
   0.038 → 0.027 → 0.020 min⁻¹ (0-30 / 30-60 / 60-120): not first order either.
4. **Ala dependence (Fig. 2F).** Excess k1(0-30) over MGO alone: 0.0043 / 0.0094 / 0.0159 / 0.0238 min⁻¹ at Ala 10 /
   25 / 50 / 100 mM — log-log slope **0.75** (0.84 on the 0-120 min window). Sublinear in Ala, assuming the MGO-alone
   channel adds. Bimolecular reading at 100 mM: k(MGO + Ala) ≈ 2.4 × 10⁻⁴ mM⁻¹ min⁻¹ (0-30) — **of the same size as
   MGO's own k2**, with Ala ~100 % recovered (catalysis, as the authors say [Q]).
5. **Temperature (Fig. 2D), k1(0-30):** Ala 1.7 × 10⁻⁴ / 4.8 × 10⁻³ / 0.044 at 60 / 80 / 100 °C → Ea(80 → 100) 121
   kJ/mol; Phe 95; Pro 97 kJ/mol. Above 100 °C MGO is gone by 30 min (no information). Apparent values: they include the
   MGO-alone channel, no come-up correction (2.5 mL ampoule in a block).
6. **Water content (SI Fig. S3B, sec. 8):** MGO left after 30 min at 100 °C is 1.2-2.7 % at 40-80 wt% water and 9 % at
   90 wt%: lower water accelerates the loss, as the authors state [Q].

## 5. What it can and cannot test in the model

**Can:**
- **The engine's MGO sink.** The shipped trunk has one MGO → melanoidin channel, `r_mgo_mel`, first order and
  amine-free, **fitted** (no literature value): log10 k(100 °C) = −1.643 ± 0.114 (0.0228 min⁻¹), Ea 20.0 kJ/mol at the
  lower search bound (`results/validation/core_prediction_uncertainty.md` rows `b1.k_mgo_mel.*`;
  `src/kinetic_core/parameters.py` `_FITTED_TEMPLATE`). Kanzler measures MGO alone at pH 5 / 100 °C as 0.010-0.014 min⁻¹
  (pseudo-first-order at 100 mM; second order in form) and with amino acids at 100 mM as 0.02-0.14 min⁻¹. **The fitted k
  sits between the amine-free and the amine-catalysed values**, which is what an amine-free first-order term fitted in
  amine pots would absorb. **The temperature dependence does not agree** (sec. 6 item 1).
- A **measured, catalytic amine dependence** of MGO loss (Ala order ~0.75-0.84, Ala not consumed), and a scavenger (Arg)
  vs catalyst (Lys, Pro) contrast.
- Lys is the only partner in common with the engine's side lanes; MGO + Lys at 100 mM: 11 % / 5.6 % left at 30 / 120
  min, Lys 98 % recovered.

**Cannot:**
- **No glycine, no cysteine**: the engine's trunk amine (Gly, `r_mgo_ak` with Gly only) and its sulfur lane (`r_mgo_cys_*`)
  get no direct test. Pro is present (the engine's `r_mgo_pro`), but only MGO loss and colour are measured, not the
  2-acetyl-1-pyrroline / acetol products.
- **pH drift unknown** (unbuffered, initial pH only); amino acids buffer weakly at pH 5 (Glu, Lys, Arg pots differ in
  buffer capacity), so partner effects and pH effects are confounded.
- **O2 unknown** (2.5 mL, ampoule headspace not stated).
- **MGO-alone order unidentified** within its SD (sec. 4.7 item 1); two time points only.
- No temperature series for MGO alone.

## 6. Conflicts and agreements with the repo's records

1. **[!] CONFLICTS with the engine's fitted `k_mgo_mel` temperature dependence.** Engine Ea 20.0 kJ/mol (clipped at the
   lower bound) vs Kanzler's apparent Ea 95-121 kJ/mol for MGO + amino acid (80 → 100 °C) and ~95-113 kJ/mol for MGO
   alone (this paper vs `bork2022`, same lab). Size at 100 °C is within the measured range; at 80 °C the engine (0.016
   min⁻¹ [Z]) is ×3.3 above MGO + Ala's 4.8 × 10⁻³ and ~×10 above MGO alone extrapolated with ~100 kJ/mol (1.6 × 10⁻³). Not registered
   as a CF (`docs/reference/` is out of scope today); proposed CF text: "MGO loss at pH 5, Kanzler 2026: apparent Ea
   95-121 kJ/mol (MGO + Ala / Phe / Pro, 80-100 C) vs fitted k_mgo_mel Ea 20 kJ/mol at its bound".
2. **[!] Form conflict (weaker):** `r_mgo_mel` is first order and amine-free; the data are second order (alone) and
   amine-catalysed (Ala order ~0.8, Ala unconsumed). A first-order term cannot carry both across concentrations.
3. **Agrees with `bork2022`** (same lab): MGO alone decelerates under a first-order reading at 130 °C too (pointwise
   k1 0.036 → 0.0087 min⁻¹), and is closer to second order there as well.
4. **[!] Internal:** the 100 °C / 30 min MGO + Ala point is 31.8 % (Fig. 1D, 2B, 2F) vs 26.5 % (Fig. 2D) vs ~60 % at
   98.4 wt% water in SI Fig. S3B (sec. 8). The first two are a ~0.2-in-ln between-set spread; the third is far outside
   and unexplained (Fig. S3 is also inconsistent in colour: 0.9 vs ~5 in Fig. 2A).

## 7. Proposed ledger entries (trunk.yml; NOT written — first reader only, second reader pending)

```yaml
  - id: kanzler2026_fig1_methylglyoxal_alone_and_with_amino_acids_ph5_100c
    dossier: kanzler2026_extraction.md
    paper: "Kanzler, Falkenhagen, Rueger, Mier & Rohn 2026 Food Chem 529:151175"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "methylglyoxal 100 mM alone or 1:1 with Ala, Val, Ile, Phe, Pro, Ser, Glu, Gln, Lys or Arg (100 mM)"
    medium: "water, no buffer, pH adjusted to 5.0 +- 0.1 with HCl/NaOH before heating; 2.5 mL in sealed glass ampoules; headspace not stated"
    temperature_c: [100]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [0, 30, 120]
    species: [methylglyoxal]
    kind: endpoint
    location: "Fig. 1D, 1E p. 4 (dossier sec. 4.1, +-1.5 % of C0)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: true
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "the pH value of the mixtures was adjusted to (5.0 +- 0.1)"
    flags: [ph_start_only, unbuffered, o2_not_stated, no_lod, below_lowest_calibrant_under_1pct, no_glycine_no_cysteine]
    note: "MGO alone 65 / 31 % at 30 / 120 min: k2 ~1.85e-4 /mM/min (second order fits). Amino acids x2.6-9.7 faster at 30 min."
  - id: kanzler2026_fig2_methylglyoxal_alanine_time_temperature_concentration
    dossier: kanzler2026_extraction.md
    paper: "Kanzler, Falkenhagen, Rueger, Mier & Rohn 2026 Food Chem 529:151175"
    lab: "TU Berlin (Kanzler / Rohn)"
    system: "methylglyoxal 100 mM with Ala 10-100 mM (also Phe, Pro 100 mM)"
    medium: "as above"
    temperature_c: [60, 80, 100, 120, 140]
    ph: {initial: [5.0], measured_during: false, final: null}
    time_min: [0, 30, 60, 120]
    species: [methylglyoxal]
    kind: time_series
    location: "Fig. 2B, 2D, 2F p. 6 (dossier sec. 4.3-4.5)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: true
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "the temperature dependence of the reaction was analyzed at a heating time of 30 min"
    flags: [ph_start_only, unbuffered, o2_not_stated, between_set_spread_100c_30min, no_mgo_alone_temperature_series]
    note: "Excess MGO loss over MGO alone ~[Ala]^0.75-0.84; apparent Ea(80->100 C) 95-121 kJ/mol vs engine k_mgo_mel Ea 20 at its bound."
```

**Provenance tiers (first reader only):** both entries `direct_measurement` (figure-digitised, external calibration),
confidence medium for the amino-acid pots (SD ±0.5-4.5 %), low for MGO alone (SD ±10 at 30 min, two time points).
Every derived k, order and Ea in sec. 4.7 is [Z], tier `derived`, not for write-back without a second reader.

## 8. Supplementary material (read 2026-11-05, first reader)

**Source.** `data/articles/kanzler2026_supplementary.docx`, copied into an empty scratch directory and read there with
`textutil` (text, tables) and a zip listing (6 EMF figures with embedded Prism 8 objects). Identified by its first lines:
this paper's title, the same five authors and address, "Supplementary material". The figures are vector EMF files;
**Fig. S3 was read from its drawing records** with a read-only EMF parser (scratchpad `tools/emf_dump.py`, maillard env,
`-I`): axis ticks and marker centres give exact plotted values [D]. Figs. S1, S2, S4-S6 hold colour, correlation and
HRMS mass-difference plots; not digitised.

### 8.1 Every figure and table

| item | caption (paraphrased) | content / values | tag |
|---|---|---|---|
| Fig. S1 | L*a*b of MGO / amino-acid pots, 100 mM, 120 min, 100 °C | colour only | — |
| Fig. S2 | log colour vs MGO conversion | the linear correlation the text states | [Q] |
| **Fig. S3** | MGO / Ala at varied water content: (A) colour, (B) MGO, (C) Ala remaining | table 8.2. [!] Caption says 120 min; the text says 30 min | [D] |
| Fig. S4 | dehydration / hydration steps per MGO in the aldol products of Table S1 | HRMS-derived | [Q] |
| Figs. S5 / S6 | frequency of mass differences in the HRMS spectra, all pots (30 min) | HRMS-derived | [Q] |
| Table S1 | ESI(−)-HRMS of MGO-only aldol products, 100 °C, pH 5, 30 min | assignments n MGO ± m H2O (e.g. 2 MGO − 0 H2O, m/z 143.0348, 39 %) | [M] |
| Table S2 | the Strecker α-aminocarbonyl of MGO + 3 MGO − 2 H2O (m/z 252.0877) in the amino-acid pots | presence only | [M] |
| Table S3 | carbohydrate-like (formose-type) species, all pots, 30 min | presence and intensity | [M] |
| Table S4 | products carrying amino-acid side chains (van Krevelen regions) | assignments | [M] |
| Tables S5 / S6 | the 25 most frequent mass differences and their repeating units | 14, 18, 30 (formaldehyde), 44 (acetaldehyde), 56 / 58 (acetol), 60 (glycolaldehyde), 72 (MGO) Da among them | [M] |

### 8.2 Fig. S3 — MGO / Ala 100 mmol/L each, 100 °C, varied water content [D, ±0.5 %; SD where drawn]

| water (wt%) | 40 | 50 | 60 | 70 | 80 | 90 | 98.4 |
|---|---|---|---|---|---|---|---|
| (A) colour A420 × F | 87.2 (7.3) | 67.8 (3.9) | 47.4 | 33.3 | 18.4 | 5.3 | 0.9 |
| (B) MGO remaining (%) | 1.6 | 1.2 | 1.3 | 1.5 | 2.7 | 9.0 | **60.2** |
| (C) Ala remaining (%) | 14.2 | 23.7 | 30.3 | 40.4 | 53.3 | 54.6 | 80.4 (10.8) |

Text: "more than 95 wt%" MGO converted at 80 wt% water; Ala conversion "around 85 wt%" at 40 wt% water; at 98.4 wt%
"a decrease of the Ala concentration could be barely observed" [M] — [!] the figure shows 80 ± 11 % left there.
[!] **The 98.4 wt% point (the paper's standard aqueous pot) does not reproduce Fig. 1D / 2B** (MGO 60 % vs 32 %, colour
0.9 vs ~5 at 30 min): either a different heating protocol for this series or a 30 / 120 min mix-up in the caption.
Use Fig. S3 only within itself (trend with water).

**Nothing in the SI on vessel volume, atmosphere, pH after heating, NF, glycine or cysteine.**

## 9. Open questions it bears on

**E16 blocker 4 (NF self-reaction with its own fragments at 20 mM) and the MGO sink.**

| item looked for | result |
|---|---|
| MGO loss alone: k, order, window | **FOUND** — 100 mM, pH 5 (initial), 100 °C, 0-120 min: k1 0.0144 → 0.0082 min⁻¹ (decelerating); **k2 ≈ 1.85 × 10⁻⁴ mM⁻¹ min⁻¹** (second order fits; order not identified within SD; two time points) |
| MGO loss with each amino acid: k, order, window | **FOUND** at 30 and 120 min for all ten (sec. 4.1, k1(0-30) 0.038-0.14 min⁻¹); time courses for Ala / Phe / Pro (not first order: Ala 0.038 → 0.020 min⁻¹); Ala order ~0.75-0.84 (10-100 mM); apparent Ea 95-121 kJ/mol (80-100 °C) |
| glycine among the ten | **NOT FOUND** |
| cysteine among the ten | **NOT FOUND** |
| lysine among the ten | **FOUND** — MGO 11 ± 4.5 % / 5.6 ± 2.5 % left at 30 / 120 min; Lys 98 ± 5 % recovered; colour ×1.7 Ala's |
| anything on NF | **NOT FOUND** — no NF pot, no NF measurement, NF not mentioned |
| stated headspace / atmosphere | **NOT FOUND** |
| pH after heating | **NOT FOUND** |

**For blocker 4 [Z]:** the paper confirms that MGO, one of NF's fragments, also fragments at pH 5 / 100 °C into
formaldehyde, acetaldehyde, acetol, glycolaldehyde and pyruvic acid (HRMS, presence only) and gives MGO's own loss
(k2 ≈ 1.85 × 10⁻⁴ mM⁻¹ min⁻¹ at 100 °C; ~2-3 × 10⁻³ at 130 °C from `bork2022`). Against `bork2022`'s NF + MGO cross
term (k2 ≳ 7 × 10⁻³ mM⁻¹ min⁻¹ at 130 °C) a free MGO fragment in an NF pot reacts with NF ≥ ~3× faster than with
itself. That supports "fragments are consumed by NF" but **does not measure NF + NF or NF + fragment in an NF-alone pot**:
blocker 4 is **not** resolved. E16 stays "undecidable on disk".

**MGO sink:** a measured constant now exists for the amine-free MGO loss at pH 5 (form second order) and for its amine
catalysis; it conflicts with the fitted `k_mgo_mel` Ea (sec. 6 item 1).
