# Buera & Karel 1993 — EXTRACTION (file stem `buera1993`)
**Buera, M. del P. & Karel, M. (1993).** "Application of the WLF equation to describe the combined effects of
moisture and temperature on nonenzymatic browning rates in food systems." *J. Food Process. Preserv.* **17**(1),
31–45. DOI **10.1111/j.1745-4549.1993.tb00224.x**.
Crossref 2026-11-03 (bibliographic query, then the DOI record): title as above, authors BUERA, KAREL, J. Food
Process. Preserv. 17(1):31-45, issued 1993-02. **Matches the PDF and the brief** (the PDF footer prints "17 (1993)
31-45"; accepted 3 March 1993). Not to be confused with **Karmas, Buera & Karel 1992** *J. Agric. Food Chem.* 40:873
(10.1021/jf00017a035, the companion paper with the Tg measurements and the Arrhenius breaks; not on disk) or with
**Buera, Chirife, Resnik & Wetzler 1987** *J. Food Sci.* 52 (the repo's `buera_1987_*` records).

### Re-analysis, no new browning rates printed. NEB rate constants from 7 literature foods (carrot, cabbage, apple, onion, potato, nonfat milk, whey powder) plus 2 own PVP model systems (glucose/glycine, xylose/lysine, 5.9 / 7.3 % moisture d.b., 45–90 °C, A420 zero order), refitted with a WLF form whose divergence temperature is set to T∞ = Tg − Tg0. Gives **Table 2: C1, C2 at To and shifted to Tg, per system and moisture**. **No k values, no Tg values, no time courses, no Arrhenius comparison are printed.**

**Read 2026-11-03 from `data/articles/buera1993.pdf`** (15 pp., 1993 scan with OCR layer). Prose from
`pdftotext -layout`. **Table 1 (p. 34) and Table 2 (p. 39) read from 220 dpi renders**, with three illegible Table 2
cells re-read from **400 dpi crops** (`ingest/buera1993_work/z_apple-09.png`, `z_ms-09.png`). Equations 2–9
(pp. 36, 37, 40, 43; not in the text layer) read from 150 dpi renders. Figs. 1–4 (pp. 38, 40, 41, 42) read from
150 dpi renders for axes, legends and form; **not digitised** (sec. 4c says why).

**Legend.** **[V]** verbatim (< 15 words) · **[M]** printed number · **[F]** read off a figure (read error stated) ·
**[Q]** authors' qualitative claim · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

**Laboratory.** Department of Food Science, Cook College, **Rutgers University**, New Brunswick NJ (not MIT; Karel had
moved from MIT to Rutgers). Corresponding author **Marcus Karel** (footnote 2); Buera a CONICET (Argentina) member on
an external grant (footnote 1). **Not in `data/keys/laboratories.yml`: NEW.** Proposed id **`karel_mit_rutgers`**
with strings `["Rutgers (Karel)", "MIT (Karel)"]`: the trunk ledger already uses the unmapped string "MIT (Karel)"
for `eichner1972_*` (2 entries). Merging the MIT and Rutgers periods of one PI is the conservative choice for
leave-one-lab-out (fewer independent labs); a split (`rutgers_karel` alone) is defensible since the methods differ
(DSC, PVP glasses vs 1972 glycerol syrups). Owner's call. Same-lab dossiers: `eichner1972_extraction.md` (MIT, Karel).
Same-lab paper cited across the repo but not on disk: Karmas, Buera & Karel 1992 (in `miao2004`, `lievonen2002`,
`k4b_paired_thresholds_and_browning.md`).

## 1. What this paper answers for the model

| model need | what Buera & Karel 1993 give | tag |
|---|---|---|
| browning rate constants vs T at low moisture | **none printed**; only log(k/k0) in Figs. 1–2 (normalised per system and moisture) | [F] |
| Tg vs moisture for foods | **measured (DSC, triplicate) but never printed**; recoverable only through the C2 identity (sec. 4b) | [Z] |
| WLF C1, C2 (fitted vs universal) | Table 2: C1⁰ 6.8–8.8 (fitted slope), C2⁰ 80–153 °C, C1ᵍ 7.2–25.1, C2ᵍ 50–120 °C; universal 17.44 / 51.6 rejected [Q] | [M] |
| reference temperature | To chosen inside each data range (30–90 °C); coefficients shifted to Tg by Eqs. 6–8 | [M] |
| T − Tg span of the data | ~0–120 K, one milk series possibly −31 K (below Tg) | [Z] |
| WLF vs Arrhenius goodness of fit | **not compared**; only "r > 0.98" for the linearised WLF plots | [M] |
| moisture effect at fixed T | **not given**: k0 is per system and moisture, so the rate level is not predicted | [Z] / [!] |
| model-system composition | Table 1: PVP 99 / glucose 0.5 / glycine 0.5 and PVP 98 / xylose 1 / lysine 1 (% dry) | [M] |
| pH, buffer, O2 | not stated (unbuffered freeze-dried solids, stoppered vials in air) | [!] |

**Answers.**
1. **Systems** [M]: literature NEB data for carrot (Legault 1951), potato (Hendel 1955), cabbage (Mizrahi 1970),
   apple (Resnik & Chirife 1979), onion (Samaniego 1991), nonfat milk (Flink 1974 and Franzen 1990), whey powder
   (Saltmarch 1981), plus two own PVP model systems at a_w 0.22 (Table 1). Moisture 0–20 % d.b. (Table 2). Tg of the
   vegetables from Karmas 1992 (DSC on freeze-dried samples equilibrated at 0–75 % RH), of dairy products from the
   lactose state diagram (Roos & Karel 1991a) [C]. Composition of the foods not given.
2. **Browning rates** [M]: model systems A420 after dilution with water (spectrophotometer range 0.02–2), linear in
   time, "zero order kinetics were assumed" [V]; literature data also zero order [Q]. **Eq. 2 labels k "(1/h)"**,
   but a zero-order A420 slope is absorbance per h [!]; units cancel in k/k0. No k, no time axis, no dilution factor.
3. **C1, C2: fitted, not universal** (confidence ~95 % that this is what the table holds). Only **C1⁰ is fitted**
   (one slope per system; mean 7.9 [M]). **C2⁰ is not fitted**: Eq. 4 gives C2⁰ = To − T∞, and T∞ is *set* to
   Tg − Tg0 with Tg0 = 50 (vegetables), 65 (potato), 100 (dairy), 90 / 120 (PVP systems). Hence **C2ᵍ = Tg0 by
   construction** and **C1ᵍ = C1⁰·C2⁰ / Tg0 is algebraic**; its "almost linear" rise with moisture (Fig. 3) follows
   from Tg falling with moisture at fixed C1⁰, not from the data [Z] (confidence ~90 %). C1ᵍ 7–25 brackets the
   universal 17.44; C2ᵍ 50–120 vs universal 51.6 (vegetables happen to sit at 50).
4. **Reference temperature**: To 30–90 °C, said to be "in the middle" of each range (Table 2 title), but apple
   (55 of 55/83) and M.S. I (45 of 45/90) sit at the lower end [!].
5. **T − Tg range**: from the derived Tg (sec. 4b) and the printed T ranges, about **0–120 K** above Tg overall;
   whey 1–22 K, PVP systems 0–50 K, vegetables 15–120 K; milk at 0 % moisture has Tg ≈ 101 °C against a 70–130 °C
   range, so part of that series may lie **below Tg**, where WLF does not apply [Z] / [!] (confidence ~60 %; per-row
   T ranges are not printed).
6. **WLF vs Arrhenius: not tested here** (confidence ~95 %). The paper reports only r > 0.98 for log(k/k0) vs
   (T − To)/(T − T∞). Over 20–40 K windows the two forms are barely distinguishable. The implied curvature is
   real but modest: Ea_eff = 2.303·R·C1⁰·C2⁰·T² / (T − T∞)² [Z] gives cabbage 1.4 % **260 → 167 kJ/mol across 30 → 52 °C**,
   cabbage 17.9 % only 110 → 94 kJ/mol. Arrhenius-equivalent Ea at To [Z]: cabbage 189 (1.4 %) → 99 (17.9 %) kJ/mol,
   milk 230 (0 %) → 140 (5.9 %), M.S. I 179, M.S. II 105 kJ/mol. The Arrhenius break itself is in Karmas 1992 [C].
7. **Usable as a (T − Tg) rate multiplier for trunk browning at low moisture? No, beyond a glass gate** (confidence
   ~85 %). (a) k0 is the per-moisture rate at To, so the equation predicts only the *temperature scaling* at each
   moisture, never the level vs moisture at fixed T; the conclusion that the moisture dependence is "implied by the
   inclusion of Tg" [V] is an overclaim (sec. 5). (b) No absolute rates to anchor anything. (c) The form saturates:
   log k → log k0 + C1⁰ as T → ∞, so at cooking T − Tg ≥ 100 K it flattens toward an Arrhenius law with an Ea far
   below the 1.4 % value; nothing here supports extrapolation from storage (20–90 °C) to frying or baking crusts
   (ref. Ferry, quoted as caution [Q]). (d) The data sets are storage of dehydrated vegetables and powders, not
   protein matrices. **What is usable**: (i) the qualitative case that a Gordon–Taylor Tg *warn/refuse gate* (B54) is
   the right instrument, and that universal WLF constants should not be used (agrees with the plan's
   "universal WLF fail by 10^3"); (ii) Table 2 coefficients as a **diagnostic prior** on the apparent Ea of
   low-moisture browning within ~0–50 K of Tg (Ea_eff 100–260 kJ/mol), if B54 ever needs one.

## 2. Conditions (own model systems; literature systems not described beyond the citations)
- **Reactants** [M] (Table 1, % of dry sample): M.S. I PVP (MW 10 000) 99, glucose 0.5, glycine 0.5;
  M.S. II PVP (MW 24 000) 98, xylose 1, lysine 1. Low concentrations chosen "to enhance diffusion effects" [Q].
  Per vial [Z]: 0.5 mL of 20 % w/w solids → ~0.1 g solids; M.S. I 2.8 µmol glucose + 6.7 µmol glycine (1 : 2.4);
  M.S. II 6.7 µmol xylose + 6.8 µmol lysine. Molality in the sorbed water [Z]: M.S. I 0.47 (Glc) / 1.13 (Gly) mol/kg;
  M.S. II 0.91 (Xyl) / 0.94 (Lys) mol/kg. Suppliers and purity not stated.
- **Preparation** [M]: frozen 24 h at −35 °C over dry ice, freeze-dried 48 h (VirTis Benchtop 3L), P2O5 desiccator
  1 week ("anhydrous"), then equilibrated over saturated K acetate, **a_w 0.22**; moisture 5.9 % (I) and 7.3 % (II)
  d.b. gravimetric.
- **Buffer / pH**: none; pH not measured [!]. **Vessel**: 2 mL glass vials, hermetically stoppered; headspace air.
- **Temperatures** [M]: 45, 55, 60, 65, 70, 80, 90 °C (incubator type not stated). **Times**: not stated [!].
- **Replicates**: DSC triplicate; browning replicates not stated [!].

## 3. Analytical methods
- **Tg** [M]: DSC Mettler TA 4000 / DSC30S, 40 µL Al pans, 5 °C/min, onset of the Cp step, second scan, triplicate
  mean. **Values not printed.**
- **Browning** [M]: A420 (HP 8452A diode array) after dilution with water into 0.02–2; zero-order slopes taken as k.
  Absolute colour; no extinction coefficient, no mass normalisation stated. ± never reported.
- **Coefficients**: Ferry's reduced-variable method: log(k/k0) = C1⁰(T − To)/(T − T∞) (Eq. 5), slope = C1⁰;
  shift to Tg by δ = To − Tg, C1ᵍ = C1⁰C2⁰/(C2⁰ − δ), C2ᵍ = C2⁰ − δ (Eqs. 6–8). No standard errors.

## 4. Data
### 4a. Table 2 (p. 39) [M], 220 dpi render; 400 dpi crops for the flagged cells
T∞ = Tg − Tg0; "T.r." = temperature range of the data used; M.C. = moisture % dry basis.

| system | To (°C) | T.r. (°C) | T∞ | M.C. | C1⁰ | C2⁰ (°C) | C1ᵍ | C2ᵍ (°C) |
|---|---|---|---|---|---|---|---|---|
| Apple | 55 | 55/83 | Tg − 50 | 1.3 / 2.2 / 5. / 8.7 / 11. / 17. | 8.8 | 83 / 103 / 112 / 118 / 129 / 143 | 14.6 / **13.? [!]** / 19.7 / 20.7 / 22.7 / 25.1 | 50 |
| Cabbage | 45 | 30/52 | Tg − 50 | 1.4 / 2.1 / 3.2 / 5.6 / 8.9 / 11.7 / 17.9 | 7.8 | 80 / 90 / 94 / 103 / 115 / 121 / 153 | 12.5 / 14.1 / 14.7 / 16.1 / 18.0 / 18.9 / 23.9 | 50 |
| Carrot | 43 | 20/50 | Tg − 50 | 5.4 / 6.2 / 8.0 | 7.4 | 98 / 103 / 108 | 14.6 / 15.3 / 16.1 | 50 |
| Milk (Flink 1974) | 90 | 70/130 | Tg − 100 | 0. / 1.2 / 5.9 | 8.1 | 89 / 125 / 146 | 7.2 / 10.1 / 11.8 | 100 |
| Milk (Franzen 1990) | 90 | 70/130 | Tg − 100 | 3. / 4. / 5. | 6.8 | 140 / 145 / 150 | 9.5 / 9.9 / 10.2 | 100 |
| Onion | 30 | 20/40 | Tg − 50 | 5.6 / 8.9 / 18.9 | 8.8 | 88 / 100 / 138 | 15.9 / 18.1 / 24.5 | 50 |
| Potato | 50 | 40/80 | Tg − 65 | 4.9 / 9.4 / 15. / 20. | 7.9 | 85 / 95 / 120 / 130 | 10.4 / 11.6 / 14.6 / 15.8 | 65 |
| Whey powder | 35 | 30/40 | Tg − 100 | 5.9 / 8.0 | 8.4 | 106 / 117 | 9.0 / 9.9 | 100 |
| M.S. I (PVP10/Glc/Gly 99:0.5:0.5) | 45 | **45/90 [!]** | Tg − 90 | 5.9 | 8.3 | 90 | 8.3 | 90 |
| M.S. II (PVP24/Xyl/Lys 98:1:1) | 55 | 45/90 | Tg − 120 | 7.3 | 6.9 | 135 | 7.8 | **120 [!]** |

Cell notes: apple 2.2 % C1ᵍ prints "13." with the last digit broken; Eq. 7 gives 8.8 × 103 / 50 = **18.1** [Z], and
Fig. 3 shows an apple point near 17–18 at ~2 % [F, ±1], so the print is most likely 18.1 (the "8" damaged).
M.S. I T.r. partly illegible ("45/9?"), consistent with the stated 45–90 °C incubations. M.S. II C2ᵍ illegible;
**120** from T∞ = Tg − 120 and from Eq. 7 (6.9 × 135 / 120 = 7.8, the printed C1ᵍ) [Z].

### 4b. Derived Tg (never printed) [Z]
From C2⁰ = To − T∞ = To − Tg + Tg0: **Tg = To + C2ᵍ − C2⁰**. Exact under the paper's own convention.

| system | Tg (°C) at the M.C. rows of 4a |
|---|---|
| Apple | 22 / 2 / −7 / −13 / −24 / −38 |
| Cabbage | 15 / 5 / 1 / −8 / −20 / −26 / −58 |
| Carrot | −5 / −10 / −15 |
| Milk (Flink) | 101 / 65 / 44 |
| Milk (Franzen) | 50 / 45 / 40 |
| Onion | −8 / −20 / −58 |
| Potato | 30 / 20 / −5 / −15 |
| Whey powder | 29 / 18 |
| M.S. I / M.S. II | 45 / 40 |

Milk 0 % → 101 °C matches anhydrous lactose (~101 °C, Roos & Karel) — the convention reproduces a known value.

### 4c. Figures (read for form, not digitised)
- **Fig. 1** (p. 38, rotated): log(k/k0) vs (T − To)/(T − T∞), T∞ = Tg − 50, for carrot (To 43), cabbage (To 45),
  apple (To 55), onion (To 30), x span −0.3 to 0.3, y −2 to 3. Legends list moisture levels that **differ from
  Table 2** (sec. 5). Cabbage 1.4 % points scatter visibly off the line.
- **Fig. 2** (p. 40): all vegetables and dairy on one line, slope ≈ 7.7 [F, ±0.5] vs the stated mean 7.9.
- **Fig. 3** (p. 41): C1ᵍ vs moisture (0–20 % d.b.) — a re-plot of Table 2.
- **Fig. 4** (p. 42): C1ᵍ vs To (0–120 °C), flat per system [Q]; the M.S. I / potato symbols are ambiguous.
- Why not digitised: every point is a ratio k/k0 with k0 per system and moisture, unrecoverable from the paper;
  the absolute rates live in the eight primary sources, which should be ledgered instead if wanted.

## 5. Flags and internal consistency
| # | flag | where |
|---|---|---|
| BK-1 | **Sign of Eqs. 2, 5, 9**: printed log(k/k0) = **−**C1⁰(T − To)/…, with C1⁰ > 0, which makes browning slower at higher T; Figs. 1–2 have positive slopes. Sign carried over from the viscosity form (k ∝ 1/η flips it). Confidence ~95 % that it is a typo. | pp. 36, 37, 43 |
| BK-2 | **Eq. 9 denominator printed "T + Tg + Tg0"**; Eq. 5 with T∞ = Tg − Tg0 requires T − Tg + Tg0. | p. 43 |
| BK-3 | **C2ᵍ = Tg0 is a definition, not a finding**; C1ᵍ's linear moisture trend is algebraic (sec. 1, answer 3). | pp. 37, 40–41 |
| BK-4 | **k0 per moisture**: the general equation does not predict moisture effects at fixed T; conclusion (g) overclaims. Confidence ~85 % (k0 definition is per T0 only, but one line through the origin for all moistures in Figs. 1–2 needs it). | pp. 36, 43 |
| BK-5 | Tg measured (triplicate DSC) but not printed; only recoverable through BK-3's identity. | pp. 33–34 |
| BK-6 | k unit "(1/h)" for a zero-order A420 slope. | p. 36 |
| BK-7 | Fig. 1 moisture legends vs Table 2: cabbage 1.4/2.1/3.2/5.9/8.9/11.8/17.4 vs table 1.4/2.1/3.2/5.6/8.9/11.7/17.9; carrot 6.4 vs 6.2; apple has no 17 % series in Fig. 1; onion legend in a_w (0.33/0.43/0.59) not moisture. | p. 38 vs p. 39 |
| BK-8 | Conclusion (d): Tg0 "between 45 and 124C"; Table 2 C2ᵍ spans 50–120. | p. 42 |
| BK-9 | "To in the middle of the range" fails for apple and M.S. I (lower end). | p. 39 |
| BK-10 | Onion C1ᵍ does not reproduce from C1⁰ = 8.8 (gives 15.5 / 17.6 / 24.3 vs printed 15.9 / 18.1 / 24.5; implies ~9.0); all other rows reproduce within ±0.1. | p. 39 [Z] |
| BK-11 | Milk (Flink) at 0 % moisture: Tg ≈ 101 °C vs T.r. 70–130 °C; part of the series may be below Tg. | sec. 4b |
| BK-12 | Fig. 1 caption cites Samaniego 1990; the reference list says 1991. Fig. 3 caption "T∞ = zt_g − Tg0" typo. | pp. 38, 41 |
| BK-13 | No ± on any coefficient, no n of points per fit, no browning replicates, no times, no pH. | throughout |

## 6. Existing records in the repo that cite this paper
`grep -rIn -i 'buera\|tb00224\|1745-4549.1993'` over `data/ docs/ src/ scripts/ tasks/` (2026-11-03): **no record
cites Buera & Karel 1993.** The DOI appears nowhere; `data/keys/papers.yml` has no entry for it. Every "Buera" hit
is a different paper, checked and set aside:
- `computational_priors.json:1176` and `benchmark_intake_registry.json:8432` `buera_1987_maillard_caramelisation_ea`;
  `papers.yml:7358`: Buera, Chirife, Resnik & Wetzler 1987 *J. Food Sci.* 52 — **not this paper**.
- `brands2001_extraction.md:180`: cites Buera 1987 (caramelisation) — not this paper.
- `miao2004_extraction.md:394, 469`: Schebor, Buera, Karel & Chirife 1999; Karmas, Buera & Karel 1992 — not this.
- `lievonen2002_extraction.md:480`, `k4b_paired_thresholds_and_browning.md:515`: Karmas, Buera & Karel 1992 — not this.
- `bell1995_extraction.md:122, 445, 450`: Buera et al. 1992; Buera, Chirife & Karel 1995 — not this.
- `papers.yml:7568` DOI 10.1111/j.1745-4549.1993.tb00730.x is Zhang, Lee & Ho 1993 (same journal and volume) — not this.
- `src/kinetic_core`: no WLF, Tg or glass-transition handling. The only a_w term is `trunk_conditions.py`
  `AW_MULTIPLIER_TABLE` (Pereyra Gonzales 2010, on `k_schiff`, floor at a_w 0.33 with a glass warning); the
  acrylamide lane carries a_w as metadata plus the B14 window. Nothing there depends on this paper.
- `tasks/data_restructure_plan.md:3133` (B54) says universal WLF fails "by 10^3"; it does not cite this paper, but
  the paper agrees that universal constants should not be used (sec. 1, answer 7).

**Mismatches: none (nothing to check).**

## 7. Proposed ledger entries (trunk.yml; NOT written)
Header rule: records where the data are and what they are, never the values. Two entries: the coefficient table
(the only printed numbers) and the normalised plots (figure_only; the primary sources are the better target).

```yaml
  - id: buera1993_tab2_wlf_coefficients_neb_low_moisture
    dossier: buera1993_extraction.md
    paper: "Buera & Karel 1993 J Food Process Preserv 17:31"
    lab: "Rutgers (Karel)"
    system: "dehydrated carrot, cabbage, apple, onion, potato, nonfat milk, whey powder (literature NEB data, 0-20 % moisture d.b.) + PVP10/glucose/glycine 99:0.5:0.5 and PVP24/xylose/lysine 98:1:1 freeze-dried model systems at a_w 0.22"
    medium: "low-moisture amorphous solids; model systems unbuffered in stoppered glass vials, air"
    temperature_c: [20, 130]
    ph: {initial: [null], measured_during: false, final: null}
    time_min: []
    species: ["WLF C1 and C2 of zero-order NEB rate constants (A420 or literature browning index), at To and shifted to Tg"]
    kind: other
    location: "Table 2 (p. 39, 220 dpi render; three cells from 400 dpi crops); dossier sec. 4a; derived Tg sec. 4b"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "Zero order kinetics were assumed and the slopes were taken as reaction constants."
    flags: [no_rate_constants_printed, tg_measured_not_printed, c2_set_by_convention_not_fitted, k0_per_moisture_no_moisture_prediction, equation_sign_typo, ph_not_stated, no_uncertainties]
    note: "Only C1 at To is fitted (mean 7.9); C2g = Tg0 by construction (T_inf = Tg - Tg0), so C1g's moisture trend is algebraic. T - Tg ~0-120 K, storage 20-130 C. Diagnostic prior on low-moisture apparent Ea (Ea_eff 100-260 kJ/mol near Tg) for B54's glass gate; not a rate multiplier (no absolute k, k0 per moisture). arms: 9 systems x 1-7 moistures (different labs' primary data, one re-analysis)."
  - id: buera1993_fig1_fig2_normalised_wlf_plots
    dossier: buera1993_extraction.md
    paper: "Buera & Karel 1993 J Food Process Preserv 17:31"
    lab: "Rutgers (Karel)"
    system: "carrot, cabbage, apple, onion (Fig. 1); + potato, nonfat milk, whey powder (Fig. 2)"
    medium: "dehydrated foods, literature storage data"
    temperature_c: [20, 130]
    ph: {initial: [null], measured_during: false, final: null}
    time_min: []
    species: ["log10(k/k0) of NEB zero-order rate constants vs (T - To)/(T - T_inf)"]
    kind: ratio
    location: "Fig. 1 (p. 38, rotated) and Fig. 2 (p. 40); 150 dpi renders, read for form only"
    status: figure_only
    quant_class: other
    simulable: false
    condition_quote: ""
    flags: [k0_per_system_and_moisture_unrecoverable, legend_moistures_differ_from_table2, replot_of_primary_sources]
    note: "Ratios only; the absolute rates are in Legault 1951, Mizrahi 1970, Resnik & Chirife 1979, Samaniego 1991, Hendel 1955, Flink 1974, Franzen 1990, Saltmarch 1981 (none on disk). Ledger those if low-moisture browning k are needed."
```
`lab` string "Rutgers (Karel)" is unmapped until the owner adds `karel_mit_rutgers` (or `rutgers_karel`) to
`laboratories.yml`.

## 8. Supplementary material
None referenced (1993 journal; no appendix, no deposited data). Nothing model-relevant is elsewhere except the
primary browning data in the eight cited sources and the Tg curves in Karmas, Buera & Karel 1992
(10.1021/jf00017a035), which would supply the Tg values this paper omits. **SI needed: no.**
