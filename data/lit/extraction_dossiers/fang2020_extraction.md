# Fang, Xiang, Sun-Waterhouse, Cui & Lin 2020 — EXTRACTION (file stem `fang2020`)
DOI: 10.1021/acs.jafc.9b06046 (printed on every page footer). **J. Agric. Food Chem. 68(6):1691-1697**; received
24 Sep 2019, published online 17 Jan 2020. Crossref 2026-11-02: "Enhancing the Usability of Pea Protein Isolate in
Food Applications through Modifying Its Structural and Sensory Properties via Deamidation by Glutaminase", authors
Luyi Fang, Huan Xiang, Dongxiao Sun-Waterhouse, Chun Cui, Junjie Lin; JAFC 68 (6) 1691-1697; issued 2020-01-17,
print 2020-02-12. **Matches the PDF.**
[!] **Enzyme naming: the brief says "protein-glutaminase"; the paper uses L-glutaminase, EC 3.5.1.2, 100 GTU/g
(Amano)** — the free-glutamine amidohydrolase, not Amano's protein-glutaminase (EC 3.5.1.44) used by
`suppavorasatit2012` and `shen2022`. See flag 1.
**Not to be confused with** `fang2009_extraction.md` / `fang2010_extraction.md` (Fang X & Schmidt-Rohr, Iowa State,
melanoidin NMR). No naming collision with the stem.

### Commercial pea protein isolate (isoelectric, 74.5 % protein, N × 5.4) dispersed 1:10 w/w in water, pH set to 9.0 (NaOH), L-glutaminase at 1:500 w/w enzyme:PPI, 50 °C shaking water bath for 3 / 6 / 9 / 12 / 24 h; then pH → 7.0, 90 °C 10 min (enzyme kill), freeze-dried. Degree of deamidation (DD, %) by Conway microdiffusion vs time (Fig. 2); SDS-PAGE, NSI vs pH, turbidity, particle size, H0, ζ-potential, FT-IR; descriptive sensory (14 panelists, 0-7 scale) of 10 % dispersions, untreated vs 12 h. **One dose, one temperature, one pH. No ammonia in absolute units, no total-amide content printed, no volatiles, no hexanal, no Maillard step, no mock-processed (no-enzyme) control, pH after incubation not reported.**

**Read 2026-11-02 from `data/articles/fang2020.pdf`** (7 pp., ACS born-digital). Prose from `pdftotext` (plain and
`-layout`); all 7 pages rendered at 150 dpi (scratchpad `ingest/fang2020_work/p150-*.png`) and read by eye; Table 1
(p. 1694) and Table 2 (p. 1696) read from those renders. **Fig. 2 (p. 1693) read from a 400 dpi crop**
(`fig2_crop.png`; axis 0-70 %, 102.4 px per 10 %; read error ±0.5 % DD). **Fig. 3A (p. 1694) read from a 300 dpi
crop** (`fig3a_crop.png`; 6.54 px per NSI unit; read error ±1 NSI %). SI not on disk.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by
me, with stated read error · **[Q]** qualitative claim by the authors · **[C]** cited (from another paper or a
textbook) · **[P]** proposed, not shown · **[Z]** derived by me, never printed · **[!]** flag.

**Laboratory.** School of Food Science and Engineering, **South China University of Technology (SCUT)**, Guangzhou;
corresponding author **Chun Cui** (also Guangdong Weiwei Biotechnology Co., Ltd.); senior co-author **Dongxiao
Sun-Waterhouse**. Method lineage: refs. 3 (Hu, Sun-Waterhouse, Liu, Cui, Wang 2019, glutaminase-treated rice protein)
and 19 (Xiang, Sun-Waterhouse, Cui et al. 2018, glutaminase-treated SPI) are the same group; ref. 12 (Liao, Zhao M,
Ren, Zhao H, **Cui C**, Hu 2010) is the DD method source and links Cui to the Zhao Mouming group.
**Not in `data/keys/laboratories.yml` (no SCUT group there at all): NEW**, proposed id **`scut_cui_c`**, lab string
"SCUT (Cui C / Sun-Waterhouse D)". Ledgers already carry two other SCUT strings, "South China University of
Technology (Yu)" (`guan2012`) and "South China University of Technology (Sun, Zhao)" (`li2020`); the registry rule
keeps groups in one university apart. Cui has co-authored with Zhao M (ref. 12), so merging with "(Sun, Zhao)" is the
conservative leave-one-out alternative (P(should merge) ~25 %). No same-lab dossier on disk.

## 1. What this paper answers for the model

| model need | what Fang 2020 gives | tag |
|---|---|---|
| ammonia released by deamidation, mmol/g protein or mM | **not printed**; only DD (%) = released / total ammonia; neither numerator nor denominator printed | [!] |
| DD vs time (enzymatic, 50 °C, pH 9.0 start) | Fig. 2: 0, 3, 6, 9, 12, 24 h; text: 17.89 % at 0 h, 59.61 % at 24 h; abstract 56.1 % at 12 h | [M] / [F] |
| DD vs enzyme dose, temperature, pH | **one dose (1:500 w/w), one T (50 °C), one starting pH (9.0)** | — |
| protein concentration in the reaction | PPI : water 1 : 10 (w/w), 74.5 ± 0.2 % protein (dry basis, N × 5.4); moisture not printed → ~69-72 g protein per L water | [M] / [Z] |
| total amide content of the isolate | measured (2 M HCl, 121 °C, 3 h, microdiffusion) **but not printed** | [!] |
| pH during / after deamidation | **not reported**; adjusted to 7.0 "with 2 M NaOH or 2 M HCl" afterwards | [!] |
| ammonia removal (dialysis, washing) | **none described**; dispersions freeze-dried directly | [M] |
| volatiles / hexanal / flavour binding | **none**; only descriptive sensory ("beany flavor", 0-7 scale) | [M] |
| heating / Maillard step | only 90 °C, 10 min, pH 7 enzyme kill; no sugar added; no browning measured | [M] |
| charge change | ζ-potential −24.13 → −32.40 (12 h) → −34.97 mV (24 h) (Table 1) | [M] |

**Answers.**

1. **DD vs dose / time / T: a time course only** (confidence ~95 %). Fig. 2 gives DD at 0, 3, 6, 9, 12, 24 h at one
   dose, 50 °C and starting pH 9.0 [F]: 18.0, 38.6, 42.8, 51.2, 56.2, 59.6 % (printed endpoints 17.89 % and 59.61 %
   [M]; abstract 56.1 % at 12 h [M]). **A rate can be read, crudely:** 20.7 DD points in the first 3 h (6.9 %/h),
   then 1.4, 2.8, 1.7 %/h over the next three 3-h steps and 0.28 %/h over 12-24 h [Z]. The 3→6 h step (4.2) is smaller
   than the 6→9 h step (8.5) — non-monotone, letters C vs B, so batch scatter is at least ~4 DD points. If forced to an
   exponential approach to a ~60 % plateau from the 17.9 % baseline, apparent k ≈ 0.15-0.22 h⁻¹ [Z] (unreliable: the
   plateau is assumed, not measured). **Ammonia is not measured as an amount** (confidence ~99 %): DD is a ratio of
   titrated ammonia (microdiffusion of the sample in saturated NaOH, 12 h) to total ammonia after complete acid
   deamidation; neither titre is printed. **The 17.89 % at 0 h is an assay baseline**, which the authors attribute
   to deamidation during the 12-h standing in saturated NaOH [Q]; it may also include ammonium already in the
   commercial isolate. The enzymatic increment is therefore DD − 17.9: 38.2 points at 12 h, 41.7 at 24 h [Z].
   **Protein in the reaction medium:** 1 g PPI per 10 g water; at 74.5 % protein (dry, N × 5.4) and an assumed 4-8 %
   moisture [P], **~69-72 g protein (N × 5.4) per L of water**, ~62-65 g/kg dispersion [Z]. **Ammonia in mM [Z/P]:**
   the paper's total-amide denominator is missing, so I borrow the pea isolate Asx + Glx of `jaeger2023` (9.286 and
   14.549 g/100 g DM at 81.22 % protein → 2.08 mmol per g protein (N × 6.25)) [C], re-base it to Fang's isolate
   (13.8 % N → 86.2 % protein on N × 6.25 → ~1.79 mmol Asx + Glx per g DM) [Z], and assume 40-60 % of Asx + Glx is
   amidated [P]: total amide-N ~0.72-1.07 mmol per g PPI DM, i.e. **~66-103 mM total amide-N in the pot** (central
   ~84 mM) [Z]. Then: **released NH3/NH4+ at 12 h ≈ 37-58 mM by total DD (central ~47), ≈ 25-39 mM by the enzymatic
   increment (central ~32)**; at 24 h ≈ 39-61 / 28-43 mM; initial rate ≈ 4.6-7.1 mM/h over 0-3 h [Z]. Uncertainty
   ±40 % from the amide fraction and moisture alone; the paper itself could have pinned it and did not.

2. **Flavour / volatiles: no volatile measured** (confidence ~99 %). The only flavour data are descriptive sensory
   (Table 2): beany flavour 5.74 ± 0.21 (untreated) vs 4.10 ± 0.32 (EPPI-12h), bitterness 1.90 ± 0.26 vs 1.20 ± 0.14,
   0-7 scale, 10 % SPI anchored at 5.0 for beany [M]. No hexanal, no GC, no binding experiment. The authors put the
   off-flavour down to beany notes plus oxidation of residual unsaturated fatty acids [Q]. **The comparison is
   confounded** (confidence ~85 %): "EPPI-0h" is the untreated commercial powder, not a mock-processed control, so
   the 12-h sample also went through pH 9 / 50 °C / 12 h, 90 °C / 10 min and freeze-drying, any of which can strip or
   transform volatiles; P(a no-enzyme processed control would show part of the beany drop) ~70 %. **No Maillard
   heating** (confidence ~95 %): the only heat step is the 90 °C, 10 min enzyme kill at pH 7 with no sugar added;
   residual sugars of the isolate not reported; no colour change on the panel (colour 4.16 vs 3.96, n.s.) [M].

3. **Could this set an initial ammonia (or amide-N) budget for a deamidated pea-protein pot? Only as an
   assumption-laden bracket, not as a measured input** (confidence ~80 %). What the engine has today (grep 2026-11-02):
   `AMMONIA` exists only in `data/species/literature_structures.yml` with the note "the engine has no ammonia
   species"; `src/kinetic_core/` has no NH3 / NH4+ species; the sulfur lane writes NH3 as a leaving product of
   cysteine Strecker / thermolysis steps and tracks no fate; the acrylamide lane's `k_asn_asp` ("Asn -> Asp + NH3",
   De Vleeschouwer 2009, free Asn, a_w 0.92) is free-amino-acid deamidation in a low-moisture model and does not
   touch protein-bound amides; `ph_state.py` `AMINE_FATE_BASIS` declares liberated amino nitrogen **non-titratable at
   release**. So there is no consumer for an ammonia budget; a deamidated-isolate pot would need (i) an `NH3` species
   with its reactions (sugar → imidazoles / pyrazine N, dienal → 2-pentylpyridine; `kim1998`, `zhou2000`) and (ii) an
   initial charge. This paper could inform (ii) as DD × amide-N × protein loading (answer 1: ~25-60 mM at ~70 g/L
   protein), with these caveats:
   - **Enzymatic vs thermal.** 50 °C, pH 9.0 start, 0.2 GTU per g PPI [Z], hours; thermal deamidation during cooking
     is a different, slower process (the repo's unverified `zhang_1993_protein_deamidation_ammonia` prior says
     0.015 h⁻¹ at 80 °C, pH 7, i.e. ~0.25 % in this paper's 10-min kill step [Z]). Nothing here gives a thermal rate.
   - **Enzyme attribution.** No no-enzyme pH 9 / 50 °C control; mild alkali at 50 °C for 24 h deamidates too, and
     EC 3.5.1.2 is a free-Gln enzyme (flag 1). The DD rise is "deamidation during the treatment", not "by glutaminase".
   - **Ammonium vs free NH3.** NH4+ pKa 9.25 at 25 °C, ΔH ≈ 52 kJ/mol [C, textbook] → pKa ≈ 8.55 at 50 °C [Z]: at
     pH 9.0 / 50 °C ~74 % is free NH3 [Z], volatile, and the vessel closure in the shaking bath is not stated, so some
     may have escaped. At pH 7 / 25 °C 0.56 % is free; at a hot cooking pot (pKa ≈ 7.0 at 120 °C [Z]) and pH 6.5 about
     22 % [Z]. pH drift during incubation is unreported (each event releases ~0.7 H+ at pH 9 [Z]).
   - **Removal before use.** No dialysis or washing: the neutralised dispersion was freeze-dried, so ammonium
     chloride (when HCl was used) stays; ammonium carboxylates can partly lose NH3 under vacuum. P(EPPI powder retains
     ≥ 50 % of the released ammonium) ~60 % [P]. If retained, a 10 % w/w dispersion of EPPI-12h (the sensory sample)
     carries ~40-65 mM NH4+ [Z] — a large base pool next to a pot's Maillard acids.
   - **Basis.** Total amide is not printed; the bracket rests on borrowed Asx/Glx and an assumed amide fraction (40-60 %).

   **Carboxylates / pI / buffering (`ph_state.py`):** each deamidation adds one side-chain carboxylate (pKa ~4-4.5 in
   proteins [C]); the enzymatic increment adds ~0.34 mmol per g PPI DM [Z], i.e. ~+38 % over the isolate's native
   Asp + Glu free carboxyls (~0.9 mmol/g by the same borrowed composition [Z]), ~25-39 mM in the pot. ζ-potential
   becomes more negative (−24.1 → −32.4 mV) [M], consistent. **But a lower pI is not visible in Fig. 3A**: the NSI
   minimum stays at pH 4-5 for all three samples [F] — the authors' "reduced pI" is from the introduction, not shown.
   `ph_state.py` today carries **no protein buffering at all** (protein-bound Lys `LYSP` is explicitly left out of
   the charge balance; no protein carboxylates), so this matters only when protein buffering enters the pH model;
   then a deamidated isolate must be a different buffer entry (more carboxylate, possibly an NH4+ pool), not the
   native `pea_isolate`. P(this changes a shipped number now) ~0 %; P(it matters if a deamidated-isolate recipe is
   ever modelled with protein buffering) ~60 %.

4. **Ledger entries:** one, in `aroma_other.yml` (the ledger whose header holds "plant-protein ... precursor-composition
   data"; the `kim1998` ammonia-release entry also sits there): the Fig. 2 DD time course as an amide-N / ammonia
   supply data set, `simulable: false`. trunk.yml does not fit (no sugar, Amadori, dicarbonyl or browning); matrix.yml
   does not fit (no binding, partition or retention measured). The sensory table is not proposed (sec. 7).

## 2. Conditions

- **Protein** [M]: commercial pea protein isolate by isoelectric precipitation (processing chart Fig. S1, SI),
  Shandong Liuliushun Foods Co., Ltd.; **74.5 ± 0.2 % protein (micro-Kjeldahl, dry basis, N × 5.4)**. Moisture, fat,
  ash, sugars not reported. (Reference SPI: Shandong Yuwang, 75.2 % dry basis, N × 5.5 — sensory reference only.)
- **Enzyme** [M]: L-glutaminase, **EC 3.5.1.2, 100 GTU/g**, Amano Enzyme (Japan), stored 4 ± 1 °C. **Enzyme:substrate
  1:500 (w/w)**, i.e. ~0.2 GTU per g PPI [Z]. Single dose.
- **Dispersion** [M]: PPI in deionised water **1:10 (mass ratio)**, magnetic stirring; **initial pH 9.0 with 3 M NaOH**.
  No buffer.
- **Incubation** [M]: shaking water bath, **50 °C, 3, 6, 9, 12, 24 h** (samples named EPPI-6h, -12h, -24h in later
  assays; DD measured at all five times). Vessel and closure not stated [!]; headspace / O2 not stated.
- **Stop** [M]: pH adjusted to **7.0 with 2 M NaOH or 2 M HCl** (the direction implies the post-incubation pH varied;
  never reported) [!]; **90 °C, 10 min** to inactivate; **freeze-dried**; stored in sealed pouches at room T. No
  dialysis, centrifugation or washing.
- **Control** [M]: "EPPI-0h" = untreated commercial PPI (not taken through pH 9 / 50 °C / 90 °C / freeze-drying) [!].
- **Sensory samples** [M]: 10 % (w/w) dispersions of EPPI-0h, EPPI-12h and reference SPI in boiled water cooled to
  25 °C; 30 mL in lidded plastic cups; served at 25 °C; three replicates; 14 trained panelists (7 M / 7 F, 20-50 y),
  three days over one week.
- **Replicates** [M]: "triplicate measurements for each of the replicate test samples obtained from separate
  experiments"; Table 1 footnote "three independent trials". Number of independent deamidation batches per time
  point not stated explicitly.

## 3. Analytical methods

- **DD** [M]: Conway microdiffusion (ref. 17, Conway & O'Malley 1942), per ref. 12 (Liao 2010). Total ammonia: 0.5 g
  sample in 5 mL 2 mol/L HCl, sealed tube, **121 °C, 3 h** ("complete deamidation"), ammonia absorbed in 20 g/L boric
  acid and titrated. Released ammonia: 1 mL "sample solution" + 3 mL saturated NaOH in the outer chamber, 3 mL boric
  acid + methyl red / bromocresol green in the centre, **12 h**, titrated with 0.01 mol/L HCl. DD = released / total.
  Whether "sample solution" is the incubation dispersion or the redissolved freeze-dried powder is not stated [!].
  The saturated-NaOH step itself deamidates (the authors' explanation of the 17.89 % at 0 h). Acid hydrolysis also
  releases some NH3 from Ser/Thr degradation (inflates the denominator slightly) [C, general]. No LOD; the ± on DD
  (error bars in Fig. 2, about ±0.5-1 point) is presumably SD of triplicates.
- **Protein** [M]: micro-Kjeldahl, N × 5.4 (PPI); Lowry for supernatants (NSI, H0).
- **NSI** [M]: 0.5 g in 10 mL water, pH 2-9 (2 M NaOH / HCl), 30 min, 10 750 g 20 min 4 °C; soluble N / total N.
- **Turbidity** [M]: 1 % w/w, pH 7 (text), A500 × ln10 / L over 20 min.
- **Particle size** [M]: Malvern Mastersizer 2000 in water (RI 1.414 / 1.33), d3,2 and d4,3 (**units not printed**,
  presumably µm) [!]. **ζ-potential** [M]: Zetasizer Nano ZS, 0.1 mg/mL in deionised water, 25 °C (**pH not stated**) [!].
- **H0** [M]: ANS fluorescence (390 / 470 nm), 0.01 M phosphate pH 7, 0.04-5 mg/mL; initial slope (arbitrary units).
- **FT-IR** [M]: Bruker Vertex 70, KBr, 400-4000 cm⁻¹, 4 cm⁻¹, 64 scans; amide I deconvolution in PeakFit 4.0.
- **Sensory** [M]: QDA after Drake 2003, 0 (absent) - 7 (extremely intense); references: 10 % SPI set at 3.0-5.0
  per attribute, 0.08 % caffeine = 5.0 for bitterness; mean ± SD; ANOVA + Tukey, p < 0.05 (Minitab 17).

## 4. Data

### 4a. Fig. 2 (p. 1693): degree of deamidation vs incubation time, 50 °C, pH 9.0 start, E/S 1:500 [F ±0.5 %]

| time (h) | DD (%) | letter | increment over previous point (%) | letter | rate (DD %/h) [Z] |
|---|---|---|---|---|---|
| 0 | 18.0 (printed **17.89** [M]) | F | 0 | D | — |
| 3 | 38.6 | E | 20.7 | A | 6.9 |
| 6 | 42.8 | D | 4.2 | C | 1.4 |
| 9 | 51.2 | C | 8.5 | B | 2.8 |
| 12 | 56.2 (abstract **56.1** [M]) | B | 5.1 | C | 1.7 |
| 24 | 59.6 (printed **59.61** [M]) | A | 3.6 (over 12 h) | C | 0.28 |

The x axis is categorical: 24 h is plotted one step after 12 h at the same spacing as the 3-h steps, and the legend
says "DD increments in every 3 h" although the last increment spans 12 h [!]. Increments read independently agree
with the DD differences within 0.2 points [Z]. Text: highest rate in the first 3 h; the 12-24 h increment not
significantly different (p > 0.05) [Q/M].

### 4b. Derived ammonia budget (not in the paper) [Z, on assumptions marked P]

Inputs: protein 74.5 % DM (N × 5.4) → N 13.8 % DM; moisture 4-8 % [P]; Asx + Glx 2.08 mmol per g protein
(N × 6.25) from `jaeger2023` (pea isolate) [C]; amide fraction of Asx + Glx 0.40-0.60 [P]; 100 g PPI per kg water.

| quantity | low | central | high | unit |
|---|---|---|---|---|
| protein (N × 5.4) per L water | 68.5 | 70 | 71.5 | g/L |
| total amide-N per g PPI DM | 0.72 | 0.90 | 1.07 | mmol/g |
| total amide-N in the pot | 66 | 84 | 103 | mM (per L water) |
| NH3 + NH4+ released, 12 h, total DD 56.1 % | 37 | 47 | 58 | mM |
| same, enzymatic increment only (38.2 points) | 25 | 32 | 39 | mM |
| NH3 + NH4+ released, 24 h, total DD 59.6 % | 39 | 50 | 61 | mM |
| same, enzymatic increment only (41.7 points) | 28 | 35 | 43 | mM |
| initial release rate, 0-3 h (20.7 points) | 4.6 | 5.8 | 7.1 | mM/h |
| free NH3 fraction at pH 9.0, 50 °C (pKa 8.55) | — | 0.74 | — | — |

### 4c. Table 1 (p. 1694): physicochemical characteristics, mean ± SD, three independent trials [M]

| sample | d3,2 (unit not printed) | d4,3 | H0 (a.u.) | ζ-potential (mV) |
|---|---|---|---|---|
| EPPI-0h | 49.82 ± 2.00 A | 101.91 ± 7.63 A | 2094.37 ± 8.59 B | −24.13 ± 1.29 A |
| EPPI-12h | 22.46 ± 1.24 B | 38.87 ± 2.62 B | 2728.80 ± 6.44 A | −32.40 ± 1.59 B |
| EPPI-24h | 26.88 ± 4.33 B | 58.93 ± 12.45 B | 2794.60 ± 3.90 A | −34.97 ± 0.31 B |

### 4d. Fig. 3A (p. 1694): NSI (%) vs pH, coarse read for the pI question only [F ±1]

| pH | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|
| EPPI-0h | 15.7 | 12.2 | 8.0 | 8.0 | 14.4 | 17.0 | 17.3 | 16.8 |
| EPPI-12h | 32.3 | 23.5 | 13.6 | 12.4 | 45.7 | 62.4 | 62.2 | 61.2 |
| EPPI-24h | 26.0 | 23.9 | 17.4 | 17.9 | 26.5 | 44.8 | 48.5 | 45.7 |

Minimum at pH 4-5 in all three; no visible downward shift of the solubility minimum [Z].

### 4e. Table 2 (p. 1696): descriptive sensory, 10 % w/w dispersions, 0-7 scale, mean ± SD [M] (* = p < 0.05 within row)

| attribute | reference (set score) | EPPI-0h | EPPI-12h |
|---|---|---|---|
| colour intensity | 10 % SPI (4.0) | 4.16 ± 0.14 | 3.96 ± 0.21 |
| **beany flavour** | 10 % SPI (5.0) | **5.74 ± 0.21*** | **4.10 ± 0.32** |
| **bitterness** | 0.08 % caffeine (5.0) | **1.90 ± 0.26*** | **1.20 ± 0.14** |
| homogeneity | 10 % SPI (3.0) | 3.51 ± 0.69 | 5.52 ± 0.56* |
| dispersibility | 10 % SPI (3.0) | 3.55 ± 0.59 | 6.54 ± 0.14* |
| suspendability | 10 % SPI (4.0) | 2.54 ± 0.00 | 6.68 ± 0.75* |
| grittiness / chalkiness | 10 % SPI (4.0) | 5.50 ± 0.85* | 1.56 ± 0.52 |
| lumpiness | 10 % SPI (4.0) | 4.85 ± 0.18* | 1.30 ± 0.54 |
| overall taste liking | 10 % SPI (4.0) | 4.14 ± 0.35 | 4.98 ± 0.26* |

### 4f. Functional data not transcribed (one line each)

Turbidity (Fig. 3B; text prints t = 0 and 20 min readings for 0 / 12 / 24 h: 4.810, 5.291, 5.630 → 1.731, 3.138,
4.241); SDS-PAGE (Fig. 1, no band change); FT-IR spectra (Fig. 4) and secondary-structure fractions (text: β-sheet
34.15 → 23.95 %, antiparallel β 19.12 → 11.66 %, β-turn 22.39 → 24.58 %, α-helix 9.90 → 11.47 %; full table is SI
Table S1); dispersion photos (Fig. S2); particle-size distributions (Fig. S3).

## 5. Flags and internal consistency

1. [!] **Enzyme identity.** EC 3.5.1.2 L-glutaminase (100 GTU/g) is a free-glutamine amidohydrolase; the
   protein-deamidating enzyme sold by Amano is protein-glutaminase, EC 3.5.1.44 (used in `suppavorasatit2012`,
   `shen2022`). The authors also say glutaminase acts on Gln **and Asn** residues [Q], which protein-glutaminase does
   not (it is Gln-specific) [C]. Either the product is mislabelled in the paper or most of the "enzymatic" DD is
   non-enzymatic. Not resolvable from the paper; P(the DD rise is mainly enzymatic) ~55 %.
2. [!] **No no-enzyme control** under pH 9 / 50 °C / 24 h; all treated-vs-untreated contrasts (DD, NSI, ζ, H0,
   sensory) are confounded with alkaline incubation, 90 °C heating and freeze-drying.
3. [!] **DD baseline 17.89 % at 0 h** is an assay artefact by the authors' own reading (12 h in saturated NaOH);
   it is assumed constant across samples, which is untested (unfolded deamidated protein may deamidate faster in
   the assay).
4. [!] **Total ammonia (amide content) measured but not printed**; neither is any titre. DD cannot be converted to
   mmol/g from the paper alone.
5. [!] **pH after incubation not reported**; "adjusted to 7.0 with NaOH or HCl" implies it moved in an unknown
   direction. Unbuffered dispersion; ~74 % of released ammonia is volatile NH3 at pH 9.0 / 50 °C [Z]; vessel closure
   unstated.
6. [!] **Fig. 2 axis is categorical** (24 h at equal spacing); legend says increments "in every 3 h" though the last
   spans 12 h. Increments are consistent with DD differences [Z].
7. [!] **Non-monotone increments** (3→6 h: 4.2; 6→9 h: 8.5, letters C vs B): batch scatter of several DD points; a
   rate constant from these six points is fragile.
8. [!] **pI claim not shown**: the introduction's "reduced isoelectric point" is not visible in Fig. 3A (minimum stays
   at pH 4-5). ζ-potential pH not stated.
9. [!] **Units of d3,2 / d4,3 not printed**; H0 arbitrary.
10. [!] Table 2: suspendability of EPPI-0h **2.54 ± 0.00** — an SD of zero across 14 panelists × 3 replicates is
    implausible (typo or rounding).
11. [!] Sensory water "boiled water (25 °C)" means cooled boiled water; the sensory samples were not heated.
12. [!] Protein factor 5.4 for PPI vs 5.5 for SPI vs 6.25 in most of the corpus: per-gram-protein comparisons
    across dossiers need re-basing (N × 6.25 equivalent: 86.2 % [Z]).
13. Consistency checks that pass: abstract 56.1 % vs Fig. 2 12-h point 56.2 [F]; text 59.61 % vs Fig. 2 59.6 [F];
    turbidity numbers in the text match Fig. 3B by eye; Table 1 d4,3 values match the text.

## 6. Existing records in the repo that cite this paper

`grep -rIn` on 2026-11-02 for `jafc.9b06046`, `9b06046`, `fang2020`, `Fang` (excluding Fang X & Schmidt-Rohr 2009 /
2011 and Fang & Cadwallader 2014), `Liuliushun`, `glutaminase`, `deamidat` across `data/`, `docs/`, `src/`,
`scripts/`, `tasks/`: **no record cites this paper.** Related records (not citing it, noted for context):
- `data/lit/computational_priors.json` → `zhang_1993_protein_deamidation_ammonia` (twice: Ea ranges, and k 0.015 h⁻¹
  at 80 °C, pH 7) — thermal protein deamidation, not this paper; not checked here.
- `data/lit/measurement_ledger/aroma_other.yml` → `kim1998_fig1_ammonia_release_amide_amino_acids_oil` and
  `shen2022_tab1_pea_protein_free_nh2_gum_conjugation` (protein-glutaminase pea isolate, Kansas State).
- `data/lit/measurement_ledger/matrix.yml` → `suppavorasatit2012_tab4_soymilk_thresholds_deamidated`.
- `data/species/literature_structures.yml` → `AMMONIA` ("the engine has no ammonia species").
- `data/lit/reaction_rules.yml` l. 394 names "glutamine or asparagine deamidation" as a free-ammonia source in a recipe.
No mismatches to report.

## 7. Proposed ledger entries (NOT written into the ledger)

**`data/lit/measurement_ledger/aroma_other.yml`** (plant-protein precursor-composition data; ammonia-release entries
live here):

```yaml
  - id: fang2020_fig2_pea_protein_deamidation_degree_timecourse
    dossier: fang2020_extraction.md
    paper: "Fang, Xiang, Sun-Waterhouse, Cui & Lin 2020 J Agric Food Chem 68:1691 (10.1021/acs.jafc.9b06046)"
    lab: "SCUT (Cui C / Sun-Waterhouse D)"
    system: "commercial isoelectric pea protein isolate (74.5 % protein dry, N x 5.4), L-glutaminase EC 3.5.1.2 at 1:500 w/w enzyme:PPI"
    medium: "PPI 1:10 w/w in deionised water, pH 9.0 with NaOH, unbuffered, shaking water bath; then pH 7.0, 90 C 10 min, freeze-dried, no dialysis"
    temperature_c: [50]
    ph: {initial: [9.0], measured_during: false, final: null}
    time_min: [0, 1440]
    species: [degree of deamidation (released / total amide ammonia), ammonia]
    kind: time_series
    location: "Fig. 2 (DD and per-interval increments at 0, 3, 6, 9, 12, 24 h; digitised in dossier sec. 4a); endpoints 17.89 % and 59.61 % in text, 56.1 % at 12 h in abstract; derived mM bracket in sec. 4b"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "incubated in a shaking water bath at 50 °C for different time intervals"
    flags: [matrix_precursor_composition, matrix_food, ph_start_only, unbuffered, condition_uncertain, no_no_enzyme_control, assay_baseline_dd_17.9pct_at_t0, total_amide_not_printed, enzyme_ec_3.5.1.2_not_protein_glutaminase, "also_relevant: trunk", "also_relevant: acrylamide"]
    note: "Only DD (%), a ratio of microdiffusion titres; neither the released nor the total ammonia is printed, so mM needs a borrowed amide-N content (dossier 4b: ~25-60 mM at ~70 g/L protein, +-40 %). One dose, one T, one start pH; no mock-processed control. Ammonia not removed before freeze-drying (retention unknown). No engine consumer: no NH3 species in src/kinetic_core."
```

**Not proposed:**
- **Table 2 sensory (beany flavour, bitterness)** — a panel score with no odorant measured, no binding, no
  partition; confounded with the processing chain (flag 2); nothing in any lane computes "beany intensity". If the
  owner wants descriptive panels ledgered as directional checks (the `brewer1995` descriptive-analysis precedent sits
  in lipid.yml), it would be a `kind: other`, `quant_class: other`, `simulable: false` entry in lipid.yml; I would not.
- **Table 1 / Fig. 3 functional data** (particle size, H0, ζ, NSI, turbidity) — no lane consumes them; ζ and NSI are
  quoted in the dossier for the pH-model note only.

## 8. Supplementary material

The paper cites Supporting Information (free at pubs.acs.org/doi/10.1021/acs.jafc.9b06046): **Table S1** (secondary
structure and aggregate fractions), **Fig. S1** (PPI processing chart), **Fig. S2** (dispersion photographs, 1 min -
2 h), **Fig. S3** (particle size distributions). **Nothing model-relevant is only in the SI** (confidence ~85 %):
the SI does not, by its description, carry the total amide / ammonia titres that sec. 1 needs. Fig. S1 might state
the isolate's extraction pH and drying route (weakly relevant to the baseline ammonium); not needed. SI not on disk.
