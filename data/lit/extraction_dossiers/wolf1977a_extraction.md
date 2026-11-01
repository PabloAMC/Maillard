# Wolf, Thompson & Reineccius 1977a — EXTRACTION
DOI: 10.1111/j.1365-2621.1977.tb08420.x (not printed in the PDF; Crossref 2026-11-01: "INITIAL LOSSES OF AVAILABLE LYSINE IN MODEL SYSTEMS", Wolf, Thompson, Reineccius, J. Food Sci. 42(6), 1540-1544, 1977-11. **Matches the brief and the lead in `wolf1977_extraction.md` sec. 7.** PDF title, authors, journal, volume and pages agree.)
### Soy protein (Promine D, Central Soya; called a "concentrate" here) + glucose + microcrystalline cellulose, slurried in 0.1 M buffer (acetate pH 5.32; phosphate pH 6.0 / 7.0; boric acid pH 8.0 / 8.68), freeze-dried, equilibrated over saturated salts to a_w 0.33-0.94 (instrument-checked), sealed retort pouches, FDNB available lysine. Three experiments: reaction order (105 °C, two protein levels), phase-three reversion by re-slurrying (105 °C; plus a 1 % glucose series at 95 °C), and a **3-factor (pH, a_w, glucose) central composite design at 80 °C / 360 min and 125 °C / 1.5 min, reduced to response-surface equations for E and k_100**. **First order in lysine; E 32.9 kcal/mol (137.7 kJ/mol) at the design centre; no per-condition k table; moisture and most glucose values measured but not printed.**

**Read 2026-11-01 from `data/articles/wolf1977a.pdf`** (5 pp., pp. 1540-1544, scanned with an OCR layer; p. 1544 is
references plus the tail of an unrelated antioxidant paper). Prose from `pdftotext`; **Tables 1-4, the response-surface
equations and Figs. 1, 3-6 read from 130 dpi page renders and 300 dpi crops** (scratchpad `wolf1977a/p-*.png`,
`eq-3.png`, `figs-4.png`, `fig6-4.png`). The OCR layer garbles the equations ("lSO*pH*G", "boo", "9S°C"); every
number below is from the render. Figure read error: Figs. 4-5 +-0.01 min⁻¹, Fig. 6 +-0.001 min⁻¹, Fig. 1 +-2 % / +-1 min.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** University of Minnesota, St Paul: Agricultural Engineering (Wolf, Thompson) and Food Science &
Nutrition (Reineccius). **Existing lab id `minnesota_reineccius`** (same lab and assay as `wolf1977_extraction.md`);
not new. Presented at IFT 1976; Minnesota Agric. Exp. Stn. paper 9634.

## 1. What this paper answers for the model

| model need | what Wolf 1977a gives | tag |
|---|---|---|
| order of available-lysine loss | **first order in lysine** in phase one: equal slopes at 7.5 % and 15 % protein (zero / second order rejected at 99 %) — "The specific reaction rates for the first order curves were 0.024 and 0.026" (p. 1542), i.e. **0.024 and 0.026 min⁻¹ at 105 °C**, a_w 0.65, pH 7.0, 4 % glucose, 18 min | [M] |
| rate constant vs composition | **no table of k per condition.** Only response-surface equations for E and k_100 in coded units (sec. 3.3) and Figs. 3-6 drawn from them | [M] / [F] |
| Ea | **32.9 kcal/mol = 137.7 kJ/mol at the design centre** (pH 7.0, a_w 0.65, glucose 2 %); Case 1 range over glucose 0.32-3.68 %: **25.9-32.4 kcal/mol (108-136 kJ/mol)** [Z from the printed equation]; each E is a TWO-temperature Arrhenius (80 and 125 °C) | [M] / [Z] |
| glucose dependence | **strongly sub-proportional at <= 105 °C**: Case 1 k(3.68 %) / k(0.32 %) = 2.3 at 80 °C, 3.8 at 100 °C, 6.6 at 125 °C for an 11.5x glucose step [Z]; the text says "5.3 times" at 80 °C and "25-fold" at 130 °C [V], which neither printed case reproduces [!] | [M] / [Z] |
| pH, a_w | act only through interactions with glucose; k rises with pH above 2 % glucose and falls with pH below it (Fig. 5); no a_w maximum near 0.65 — "a maximum k at water activities near 0.65 did not occur" | [Q] / [F] |
| plateau ("phase three") | at 105 °C, 4 % glucose, available lysine 9.4 -> 6.5 (30 min) -> 4.0 (42) -> 4.0 (54) mg/g; re-slurried at 30 min it reaches 3.6 (42) and 4.9 (54), n.s. from continuous (Table 2). **Glucose is NOT exhausted**: "the glucose levels decreased from 4% to 2.5%" in the first 30 min | [M] |
| plateau at low sugar | 1 % glucose, 95 °C: 7.2 -> 6.7 (15) -> 4.9 (30) -> 5.4 (45) -> 5.0 (60) -> 4.9 (75) mg/g (Table 3); glucose not printed for this series [!] | [M] |

**Answers.**
1. **The constants the companion paper pointed to are here, but only as regression surfaces** (confidence 0.95): k_100
   = 0.0124 + 0.0043·G_coded min⁻¹ and E = 32.9 + ... kcal/mol (Case 1). At 125 °C, 4 % glucose (G_coded +2, just
   outside the +-1.68 design) Case 1 gives **k = 0.30 min⁻¹** (Case 2: 0.57) [Z], consistent with the companion's
   first-pass bound >= 0.33 min⁻¹ (`wolf1977_extraction.md` sec. 1) within 0.04 dex.
2. **Against the engine's glycation step (F05)** [Z, moisture NOT printed; scanned w = 0.08-0.12 g water / g dry]:
   the engine's k_glyc·[Glc]_water at 4 % glucose is **0.26-0.39 min⁻¹ at 125 °C (Wolf Case 1: 0.30, +-0.1 dex)**,
   but **0.036-0.054 at 100 °C (Wolf 0.021, +0.24 to +0.42 dex)** and **0.0062-0.0093 at 80 °C (Wolf 0.0019, +0.5 to
   +0.7 dex)**. The disagreement is the barrier (engine 96.8 kJ/mol, Wolf 132-138 at 2-4 % glucose) and the glucose
   order (engine 1, Wolf ~0.35-0.8 depending on T) — both independent of CF-072.
3. **CF-072 (does the Amadori decay return lysine?)**: **not decided** (P 0.9 that this paper cannot decide it). The
   only observable on which the two rivals differ by > 0.3 dex is the 105 °C plateau (Table 2), and it separates the
   lysine-returning structure only from a BARE no-return structure; a no-return structure with an inaccessible lysine
   fraction (which the engine already carries as `available_fraction(amine_band)`) fits as well. Numbers and the
   proposed decider are in the findings note (scratchpad `wolf1977a_findings.md`), summarised in sec. 5.

## 2. Conditions

- **System** [M]: soy protein (Promine D, Central Soya) — "soybean protein concentrate" (p. 1540) [V]; the companion
  calls the same product an isolate [!] — + glucose + Avicel microcrystalline cellulose + water (Table 1).
  Protein **15 % dry wt** in experiments B and C; **7.5 or 15 %** in the order experiment A.
- **Buffer / pH** [M]: dry blend slurried in **0.1 M buffer, sample : buffer 1 : 1.5 w/v** — "pH 5.32-acetate buffer,
  pH 6.0 and 7.0-phosphate buffer, pH 8.0 and 8.68-boric acid buffer" [V]. **[!] pH is confounded with buffer species**
  (acetate / phosphate / borate; borate complexes sugar diols). Buffer salt ~0.15 mmol / g dry [Z], i.e. ~1-2 M in the
  water phase at a_w 0.65 [Z, moisture assumed]. In-situ pH not measured; pH is the slurry set point.
- **Water** [M]: freeze-dried (blast-frozen -30 °C), crushed, equilibrated 4 d over saturated salts under partial
  vacuum — "The samples were allowed to equilibrate for 4 days." [V]; a_w checked with a Beckman SMT electric
  hygrometer. Design a_w **0.33 / 0.50 / 0.65 / 0.80 / 0.94** (Table 1). **Moisture measured (GC, Hollis / Reineccius-
  Addis, wet basis) but no value printed** [!].
- **Glucose** [M]: design levels **0.32 / 1.00 / 2.00 / 3.00 / 3.68 %**; 4 % in A and B; 1 % in the Table 3 series.
  Basis not stated (protein is "% dry wt"; I read glucose as dry basis too) [!]. Glucose by hexokinase / G6PDH
  (Calbiochem Glucose Stat Pack, A340); **the only glucose value printed is 4 % -> 2.5 % after 30 min at 105 °C**.
- **Heating** [M]: sealed retortable foil pouches, <= 0.07 g sample per cm² — "no more than 0.07g of sample per sq cm
  of surface area" [V]; water bath or miniature retort; cooled to room temperature. Above 100 °C, lags corrected by
  Gondo's iterative method — "the heating and cooling lags were corrected for by an iterative procedure" [V]; the
  corrected equivalent times are not printed.
- **Oxygen**: sealed pouches, headspace not stated (the slurry was not deaerated).
- **Times** [M]: order experiment 18 min at 105 °C (5 samples + 5 controls per protein level); reversion 30 / 42 / 54
  min at 105 °C, re-slurry after 30 min then 12 or 24 min more (Fig. 2, 3 samples per box); Table 3 0-75 min at 95 °C;
  design "The two processing temperatures were 125°C and 80°C." [V] with "intervals of these two temperatures were
  1.5 min and 360 min" [V].
- **Analysis** [M]: available lysine by FDNB (Carpenter 1960, Booth 1971 modification), "a subtractive procedure";
  "This method according to Carpenter measures biologically available lysine" [V]. No detection limit or recovery
  printed. **C0 for the design came from center-point controls only** — "24 control samples (prepared as center points
  but not processed)" [V] [!]: composition-specific C0 was never measured.
- **Statistics** [M]: "The 96 samples required for this design and 24 control samples" [V], random order; block
  (analysis-period) effect checked, none; outliers removed by range inspection — "Five of the 96 data points were
  discarded using this method." [V]; per-composition E and A by least squares on ln[ln(C/C0)/(-t)] vs 1/T; stepwise
  response surfaces (SPSS), terms kept if ΔR² >= 4 % and partial F >= 2.7.

## 3. Quantitative results

### 3.1 Table 1 (design levels, p. 1540) [M]
| | A (order) | B (reslurry) | C +1.68 | C +1.00 | C 0 | C -1.00 | C -1.68 |
|---|---|---|---|---|---|---|---|
| pH | 7.00 | 7.00 | 8.68 | 8.00 | 7.00 | 6.00 | 5.32 |
| a_w | 0.65 | 0.65 | 0.94 | 0.80 | 0.65 | 0.50 | 0.33 |
| glucose (%) | 4 | 4 | 3.68 | 3.00 | 2.00 | 1.00 | 0.32 |
| protein (% dry wt) | 7.5, 15 | 15 | 15 | 15 | 15 | 15 | 15 |

Coding [Z from Table 1]: G_c = (G% - 2) / 1; pH_c = pH - 7; a_w,c = (a_w - 0.65) / 0.15 (a_w at -1.68 is 0.33, at
+1.68 0.94: the axial a_w steps are 0.32 / 0.29, not exactly 1.68 x 0.15 [!]).

### 3.2 Reaction order (p. 1541-1542) [M]
- 18 min at 105 °C: loss **2.3 mg/g (7.5 % protein) and 4.1 mg/g (15 %)** — "available lysine decreased by 2.3 mg and
  4.1 mg/g sample, respectively" [V].
- First-order slopes **0.024 and 0.026 min⁻¹**; equality not rejected even at 50 %; zero / second order rejected at
  99 %. C0 of these samples not printed; 0.026 min⁻¹ over 18 min implies C0 ~11 mg/g at 15 % protein [Z], above
  Table 2's 9.4 and Table 3's 7.2 [!].

### 3.3 Response surfaces (p. 1542) and Table 4 (p. 1543) [M]
Coded variables G, pH, a_w; E printed with unit "(Kcal/mol·K)" [!: read as kcal/mol; the methods give "Cal/mole"].
- **Case 1** (linear, square, interaction): E = 32.9 + 1.94·G - 1.64·a_w·G - 1.32·G² + 1.50·pH·G;
  **k_100 = 0.0124 + 0.0043·G (min⁻¹)**.
- **Case 2** (+ 3rd / 4th order): E = 32.9 + 1.10·G³ - 1.77·a_w·G - 1.32·G² + 1.50·pH·G;
  **k_100 = 0.0122 + 0.00175·G³ + 0.00104·G·a_w (min⁻¹)**.
- **Case 3** (Case 2 with coded values shifted by +2): E = 18.32 + 15.5·G - 4.60·G² + 0.121·G⁴ + 1.50·pH·G -
  2.59·pH - 1.80·a_w·G **"= 3.3*a_w"** (sic; read as "- 3.3·a_w") [!]; k_100 = 0.00176 + 0.004320·G + 0.000890·a_w.
  **Case 3 does not reproduce Cases 1-2**: its centre E evaluates to 19.9 kcal/mol, not 32.9 [Z] — a misprint or an
  unstated shift convention. **Do not use Case 3.**
- Then k = k_100·exp[-(E/R)(1/T - 1/373.15)].
- **Table 4** [M] (overall F; R²; SD): Case 1 E 6.8, 0.73, 2.0 kcal/mol; k_100 23.8, 0.64, 0.003 min⁻¹. Case 2 E 9.6,
  0.79, 1.7; k_100 15.4, 0.72, 0.0031. Case 3 E 5.1, 0.84, 1.9; k_100 13.1, 0.69, 0.0032. "between 0.64 and 0.72 for
  k100 and 0.73 and 0.84 for E" [V]. **Residual SD of k_100 ~0.003 min⁻¹ = +-24 % (0.1 dex) at the centre** [Z].
- **Which case the figures use** [Z]: Case 2 reproduces Figs. 4-6 within read error (Fig. 4: 0.114 / 0.186 vs read
  0.11 / 0.167; Fig. 5: 0.236 / 0.486 and 0.267 / 0.130 vs 0.21 / 0.44 and 0.28 / 0.135; Fig. 6: 0.0127 / 0.0176 and
  0.0117 / 0.0068 vs 0.012 / 0.021 and 0.014 / 0.0078). Case 1 has no a_w term in k_100 and gives a flat Fig. 6.
- **k at named conditions** [Z, Case 1, pH 7, a_w 0.65]: 4 % glucose: 0.0019 (80 °C), 0.012 (95), 0.021 (100), 0.037
  (105), 0.30 (125), 0.50 (130) min⁻¹. 1 % glucose: 0.00084 / 0.0047 / 0.0081 / 0.014 / 0.10 / 0.16. Case 2 is within
  0.1 dex at <= 105 °C and 0.28 dex higher at 125 °C, 4 %.
- **Text claims not reproduced** [!]: "increases the rate only 5.3 times" at 80 °C and "the same increase in glucose
  causes a 25-fold increase" at 130 °C [V]; Case 1 gives 2.3 / 7.3, Case 2 2.4 / 15 [Z] (pH 7, a_w 0.65; other pH or
  a_w do not reach 25 within the design).

### 3.4 Figures [F]
- **Fig. 1** (schematic, "isothermal process at 110 °C", % remaining): 100 -> minimum ~47 % at ~13 min, step up to ~52
  % by ~15 min, flat to ~125 min. A drawn curve, no data points. Do not key.
- **Fig. 3** (contours of k, min⁻¹, glucose 0-4 % vs T 60-140 °C; pH, a_w not stated): contours .0001 to 1; at >= 0.8 %
  glucose the contours are nearly vertical (T dominates), below 0.8 % they bend sharply (glucose dominates).
- **Fig. 4** (120 °C, a_w 0.65, glucose > 2 %): k 0.11 (pH 5) -> 0.167 (pH 9), linear.
- **Fig. 5** (130 °C, a_w 0.87): 3.2 % glucose 0.21 (pH 5) -> 0.44 (pH 9); 0.8 % glucose 0.28 -> 0.135; crossing at
  pH ~5.7, k ~0.25.
- **Fig. 6** (100 °C, pH 7.0): 3.2 % glucose 0.012 (a_w 0.35) -> 0.021 (0.93); 0.8 % glucose 0.014 -> 0.0078; crossing
  at a_w ~0.43. The y axis has a break between 0 and -0.09.

### 3.5 Table 2 (reversion, 105 °C, 4 % glucose, a_w 0.65, pH 7.0, 15 % protein; p. 1542) [M]
Available lysine, mg remaining / g sample, mean +- SD (n = 3 per box in Fig. 2 [!: n not printed in the table]):
0 min 9.4 +- 0.7; continuous 30 min 6.5 +- 1.0, 42 min 4.0 +- 0.2, 54 min 4.0 +- 0.2; re-slurried at 30 min then
heated to 42 min 3.6 +- 0.3*, to 54 min 4.9 +- 0.6* (* "Not significantly different from the 42 and 54 min
continuous process"). The post-reslurry 30-min analysis drawn in Fig. 2 is not printed [!].
Fractions [Z]: 0.69 / 0.43 / 0.43 continuous; 0.38 / 0.52 reslurry. The 30 -> 42 min drop (-0.21 dex in 12 min) is
steeper than 0 -> 30 (-0.16 dex in 30 min): not a single exponential; read as batch noise (SD 1.0 at 30 min).

### 3.6 Table 3 (1 % glucose, a_w 0.65, pH 7.0, 95 °C; p. 1542) [M]
"glucose levels were set at 1%, water activity at 0.65 and pH at 7.0" [V]; "processing at 95°C and decreases in
available lysine levels occurred" [V] (the OCR layer reads "9S°C"; the render reads 95 °C). Available lysine, mg/g,
mean +- SD: 0 7.2 +- 1.6; 15 6.7 +- 1.6; 30 4.9 +- 2.0; 45 5.4 +- 1.3; 60 5.0 +- 0.8; 75 4.9 +- 0.2. Fractions [Z]:
0.93 / 0.68 / 0.75 / 0.69 / 0.68. Plateau ~0.69 from 30 min on. C0 7.2 vs Table 2's 9.4 at the same 15 % protein
(0.12 dex batch-to-batch) [!].

### 3.7 Cited background [C]
- Thompson et al. 1976 (Trans. ASAE, not on disk): three phases; phase one "a rapid loss where 30-60% of the
  available lysine is made unavailable" [V]; "Phase two shows a statistically significant (98% confidence) increase in
  available lysine" [V]; maximum loss "independent of temperature between 80°C and 130°C" [V] but dependent on glucose.
- Jokinen et al. 1976 (J. Food Sci. 41:816, not on disk): phase-three equation; "water activity explains only a small
  portion (4%) of the variation" [V].
- Warmbier et al. 1976 agrees on first order at lower T; Dworschak & Hegedus 1974 reported fourth order.

## 4. Mechanism evidence and its kind

Within-study factorial (pH x a_w x glucose x T) on ONE marker (FDNB lysine), reduced to regression surfaces; one
re-slurry intervention; one low-glucose time course. Kind: kinetic (apparent), not mechanistic: no Amadori /
furosine, no CML, no 3-DG, no colour. Hypotheses addressed [Q]: sugar exhaustion as the cause of phase three (rejected
at 4 % glucose: 2.5 % remains), re-mixing reverting phase three (rejected, Table 2). Not addressed: protected lysine
fraction, FDNB response to Amadori lysine, lysine return from Amadori decay. pH effects are buffer-confounded.

## 5. Bearing on open questions

- **F05 unknown "pea / soy isolate glycation KINETICS at cooking temperature"** — FILLED in part: a two-temperature
  (80 / 125 °C) soy E and k_100 surface over pH 5.3-8.7, a_w 0.33-0.94, glucose 0.3-3.7 %, a_w instrument-checked.
  Limits: endpoint k (one time per T), C0 from center controls, moisture not printed, buffer-confounded pH.
- **F05 engine mismatch (d), no a_w term on glycation** — PARTLY informs: a_w acts only through a_w x glucose
  (Case 2 k_100 term 0.00104·G·a_w; Fig. 6), and E carries -1.64 to -1.80 kcal/mol per coded a_w x G. Small (Jokinen:
  4 % of variance). Supports a weak, sugar-dependent a_w term, not a bell.
- **Engine barrier for protein glycation** — CONTRADICTS at the 0.3-0.7 dex level below 105 °C: Ea 96.8 (engine,
  Martins glycine) vs 132-138 kJ/mol (Wolf, soy, 2-4 % glucose). Matches magnitude at 125 °C.
- **CF-072** — LEAVES OPEN. Rival predictions [Z] (shipped B20 constants; and each structure refitted to Table 2):
  - Table 2 plateau, 54 min: shipped return structure 0.32-0.42 (obs 0.43: -0.13 to -0.01 dex); shipped no-return
    0.05-0.14 (-0.93 to -0.50 dex). With Wolf's own phase-one k (0.025 min⁻¹) and the shipped k_flp_decay, the
    return structure's steady free fraction is 0.46 (+0.03 dex). BUT a no-return structure refitted to Table 2
    (k_glyc 0.3-0.5 dex below shipped) fits with rms 0.06 dex, and the return structure refitted collapses to
    k_flp_decay -> 0 (identical fit): Table 2 alone has no structural power.
  - Table 3 (1 %, 95 °C), 75 min: return 0.72-0.79, no-return 0.62-0.71 (shipped), obs 0.68: rivals differ <= 0.07 dex.
  - Design end points (k_app = -ln(C/C0)/t): 125 °C / 1.5 min rivals differ <= 0.04 dex; 80 °C / 360 min 0.2-0.3 dex,
    in opposite directions at low and high glucose; both miss Wolf's glucose slope at 80 °C by ~0.7 dex.
  - Plateau depth vs glucose (0.43 at 4 %, 0.69 at 1 %) is the return structure's signature (free fraction
    k_d / (k_d + k_glyc·G)) and matches it within 0.06 dex; a fixed inaccessible pool predicts equal depths. But the
    1 % series runs 10 °C cooler, its glucose is not printed (exhaustion at ~1.1 mol glucose / mol lysine cannot be
    excluded, since at 4 % glucose fell ~4 mol per mol lysine lost [Z]), and C0 differs by 0.12 dex. Weight low.
  - Likelihood ratio, return : bare no-return ~2; return : no-return + inaccessible pool ~1. **What would decide**: a
    sugar-chase on plateau samples (heat after removing free glucose: return predicts available AND LC-MS total-free
    lysine rising at ~k_flp_decay, i.e. >= +0.2 dex in 60 min at 105 °C; no-return predicts none), or Nguyen 2016
    Fig. S1B's lysine balance (CF-072 remedy). Thompson 1976's phase-two rise (98 %) is the only in-lab hint of a
    return; it is not on disk.
- **E5 buffer catalysis** — CAUTION: pH arms change buffer species (acetate, phosphate, borate) at ~1-2 M water-phase
  buffer [Z]. Fig. 4-5 pH slopes must not be read as a pure pH law.

## 6. Disagreements found

1. **Engine barrier and glucose order vs Wolf (NEW)**: Ea 96.8 vs 132-138 kJ/mol; glucose order 1 vs Wolf's k
   ratio 2.3-6.6 for an 11.5x glucose step (80-125 °C). Engine k_glyc·G over-predicts by +0.24 to +0.7 dex at 80-100 °C,
   agrees at 125 °C. Candidate causes: plateau-capped 80 °C / 360 min end points inflate E (P ~0.5); C0 from center
   controls only; moisture-dependent water-phase glucose; real barrier difference between free glycine and
   protein-bound lysine in a solid matrix.
2. **Companion Table 1 row 5 is cited to this paper but is not in it** [!]: `wolf1977_extraction.md` sec. 3 (and the
   ledger note) list "5 once with 4 % glucose 5.0 +- 0.7 (from Wolf 1977a)" at 125 °C / 2 min. Wolf 1977a prints no
   125 °C, 2-min, 4 % glucose value (its tables are 105 °C and 95 °C; the design raw data are not printed). The
   companion's citation is unverifiable here (perhaps the 1978 paper or unpublished).
3. **Product naming** [!]: Promine D is "soybean protein concentrate" here, a "soy protein isolate" in the companion
   (and 20 % vs 15 % protein). Matters for residual carbohydrate (possible extra reducing sugar).
4. **Lysine per protein** [!]: C0 9.4 (Table 2) and 7.2 (Table 3) mg/g at 15 % protein = 6.3 and 4.8 g / 100 g
   protein; the order experiment implies ~11 mg/g [Z]; the companion 9.6-10.2 mg/g at 20 % protein (4.8-5.1 g / 100 g).
   Batch spread ~0.15 dex within one lab and one assay.
5. **Text vs equations** [!]: the 5.3x / 25x glucose effects (p. 1543) are not reproduced by Case 1 or 2; Case 3's
   centre E is 19.9, not 32.9 kcal/mol; the Case 3 E equation prints "=" for a sign.

## 7. Ledger-ready measurements

- Table 2 (105 °C, 4 % glucose, continuous and re-slurried) and Table 3 (95 °C, 1 % glucose): time series, printed
  mean +- SD, FDNB available lysine; simulable only with a declared moisture band (not printed).
- The order experiment (0.024 / 0.026 min⁻¹, 105 °C, two protein levels): a rate_constant row, endpoint-derived.
- The Case 1 and Case 2 response surfaces (E, k_100): a rate_constant row, `kind: rate_constant`, fitted
  regression, NOT a fit target (endpoint k, center-point C0, buffer-confounded pH); usable as a held-out TEST and a
  band for the protein-glycation barrier in a low-moisture plant matrix. Do NOT key Case 3 or Fig. 1.
- Proposed YAML snippets: scratchpad `wolf1977a_findings.md` sec. 4.
