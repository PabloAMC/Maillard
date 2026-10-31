# Hayashi & Namiki 1980 — EXTRACTION (the early C2 fragment of sugar + amine: glyoxal dialkylimine from D-glucose + alkylamines in refluxing ETHANOL; isolated as glyoxal dicyclohexylimine, 12 / 74 mmol per mol glucose at 1 : 1 / 1 : 5 glucose : amine; time course in arbitrary units peaking at 20 min, BEFORE 3-deoxyglucosone, radical and browning; route from the glycosylamine proposed, not shown; no water, no pH, no Amadori measurement)

### A QUALITATIVE-ORDER PAPER WITH TWO ABSOLUTE YIELDS. It establishes that an amine-bound C2 (glyoxal oxidation level) appears within minutes, ahead of 3-DG, and that it needs the amine (alkali alone gives none). The only absolute numbers are two isolated yields in ethanol; the kinetics are TLC-densitometer peak areas. Nothing is in water or at a controlled pH.

**Source on disk:** `data/articles/hayashi1980.pdf` (6 pp., Agric. Biol. Chem. 44(11) 2575-2580; e.Typist OCR of a
scan). Prose from the `pdftotext -layout` layer (word-spacing garbled but legible); **Figs. 5-7 read from 130 / 300 dpi
renders** (scratchpad `pdf3/hay-*.png`, `hay5.png`, `hay6.png`).

**Legend.** **[V]** verbatim (short quote) · **[M]** a number printed in this paper · **[F]** read from a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed / hypothesised, not
measured · **[Z]** derived by me, never printed.

## 1. What this paper answers for the model

The question (`docs/guides/MECHANISM_REVIEWS_2026-10-31.md` sec. 3, Tier 1; `docs/guides/STATUS_AND_PLAN.md` item
"Oxygen as a reagent ... the Namiki route"): is there a pre-Amadori C2 (glyoxal-level) channel from sugar + amine,
with yield, time, temperature, pH and amine dependence? The engine's only glyoxal source is oxidative glucosone
cleavage (`r_g_go`, Kocadagli 2016, amine-free) — no amine-dependent C2 route.

| what the model needs | what Hayashi & Namiki 1980 give | tag |
|---|---|---|
| identity of the C2 fragment | **glyoxal dialkylimine** (R-N=CH-CH=N-R); isolated and identified for R = cyclohexyl (mp 151-152 °C, M⁺ 220, ¹H-NMR δ 7.82 2H, IR vs authentic); for t-butyl only via its hydrolysis product glyoxal (TLC, GLC of reduced / acetylated product = ethylene glycol diacetate, bis-phenylhydrazone) | [M] |
| yield | **12 mmol / mol glucose** (glucose + cyclohexylamine 1 : 1, 1 M each, 95 % ethanol, reflux 15 min) and **74 mmol / mol** (1 : 5, glucose 1 M, 99 % ethanol, reflux 30 min) — isolated crystals, so lower bounds | [M] |
| time course | glucose + t-butylamine 1 M each, 95 % ethanol, reflux: C2 imine (peak area, arbitrary) ~4.5 at 10 min, **max ~6.7 at 20 min**, ~4.4 at 40, ~1.9 at 80, ~1.3 at 120 min | [F] |
| order vs other intermediates | "formed at the very early stage ... maximum after 20 min ... This precedes either the initial formation of 3-deoxyglucosone, that of the free radical species or the appearance of the brown color" | [V/Q] |
| 3-DG in the same pot | ~0.07 mM at 15 min, ~1.0 mM at 35 min, ~3.6 mM at ~88 min (≈ 0.36 % of the glucose by 88 min) | [F] / [Z] |
| temperature | reflux of 95-99 % ethanol only (~78 °C, not printed) | [Z] |
| pH | **not measured** (ethanol; amine as base). NaOH 0 → 1 M (with 2 M t-BuNH₂) slightly LOWERS the C2 imine (peak area ~3.2 → ~1.9) while browning stays flat (~1.1-1.2); glucose in 0.1 M NaOH alone gives no glycolaldehyde, glyoxal or glycolic acid | [F] / [Q] |
| amine dependence | rises with amine excess, **maximum at ~10-fold amine : glucose** (Fig. 6A: peak area ×1/5 ~2 at 0.5, ~10.5 at 4, ~19 at 7, ~20 at 16, ~4 at 32), browning peaks at ~2 and vanishes above ~8 | [F] / [Q] |
| amines / sugars covered | TLC (Fig. 1): glucose, fructose, arabinose + t-butylamine; glucose + n-butylamine, cyclohexylamine; ethylamine mentioned. **No amino acid, no protein amine** | [M] |
| before the Amadori compound? | **not measured.** The Amadori compound is never analysed; precedence over 3-DG is shown; formation from glucose or the glycosylamine (Route B, Fig. 7) by retro-aldol at C2-C3 is PROPOSED, favoured over the free sugar (Route A) because alkali lowers it and the amine raises it; "The evidence, however, is not yet conclusive" | [P] / [V] |

**Answers.**

1. **Quantified? Partly** (confidence ~95 %). Two isolated yields (1.2 % and 7.4 % mol/mol glucose) for cyclohexylamine
   in ethanol; the time course and amine / alkali series are relative only (densitometer peak area of the glyoxal
   spot after on-plate hydrolysis). No rate constant can be fitted; at best a share: in the 1 : 1 run ≥ 1.2 % of the
   glucose reaches the C2 diimine within 15 min of reflux [Z].
2. **Before Amadori? Probably, not shown** (P ~0.6 that it is pre-Amadori, my weighting of the authors' circumstantial
   case). It peaks while 3-DG is still < 0.1 mM, but 3-DG is not the Amadori compound and the Amadori compound was
   not measured.
3. **Transfer to the engine: poor.** Anhydrous ethanol with alkylamines (pKa ~10.5, unprotonated) is far from aqueous
   amino-acid pots at pH 5-8 (the neutral-amine fraction of an amino acid at pH 6 is ~10⁻³-10⁻⁴ [Z, from a typical
   pKa 9-10, not printed here]).
   Use as **a qualitative signature only**: an amine-dependent glyoxal-level C2 that appears early and decays, absent
   without amine. Its oxidation state: Route B draws CH=NR / CH₂OH → CH=NR / CHO with no oxidant shown — the C2 diimine
   is at the glyoxal level, so an oxidation (or a disproportionation) is implied but not discussed [Z].
4. **No supplementary data** (sec. 6).

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Two-Carbon Sugar Fragment at an Early Stage of the Browning Reaction of Sugar with Amine" |
| Authors | Tateki Hayashi, Mitsuo Namiki — Dept. of Food Science and Technology, Faculty of Agriculture, Nagoya University |
| Venue | *Agric. Biol. Chem.* 1980, **44**(11), 2575-2580; received March 26, 1980 |
| DOI | 10.1271/bbb1961.44.2575 (not printed; Crossref returns this title, these two authors, Agric. Biol. Chem. 44(11) 2575-2580, 1980; checked 2026-10-31). **Matches the brief.** |
| Paper type | Isolation / identification (TLC, GLC, IR, ¹H-NMR, MS) + one semi-quantitative time course (TLC densitometry, ESR, A420) |
| Related | Namiki & Hayashi 1975 JAFC 23:487 and Hayashi, Ohta & Namiki 1977 JAFC 25:1282 (pyrazine cation radicals, ref 1); the promised follow-up on the pathway is not on disk |

## 2. Methods as they matter to a model

- **Pots** [M]: sugar + amine 1 M each (Fig. 1, 5), 95 % ethanol, reflux 10-40 min; amine series with glucose 1 M
  (Fig. 6A, 40 min); NaOH series with glucose 1 M + t-BuNH₂ 2 M (Fig. 6B, 40 min). Preparative: glucose 0.1 mol +
  t-BuNH₂ 0.1 mol in 100 mL 95 % ethanol, 40 min; glucose 0.05 mol + cyclohexylamine 0.05 mol (50 mL 95 % EtOH, 15
  min) or 0.25 mol (50 mL 99 % EtOH, 30 min).
- **C2 assay** [M]: 1 µL of the mixture on silica TLC, 2,4-DNP spray, glyoxal spot (Rf ~0.7) read by densitometer
  at 400 nm after 10 min; the diimine hydrolyses to glyoxal on silica / acid (Fig. 3), so the assay reads the imine
  as glyoxal. **Any free glyoxal would read the same** [Z].
- **3-DG** [M]: 20× diluted mixture, NaBH₄ reduction, acetylation, GLC (XE-60, 180 °C), against authentic 3-DG.
- **Radical**: ESR intensity (arbitrary); **browning**: A420 × dilution.

## 3. Numbers printed

- Isolated glyoxal dicyclohexylimine: 0.14 g from 0.05 mol glucose (1 : 1, 15 min) → 0.64 mmol = 12.7 mmol/mol [Z,
  matches the printed 12]; 0.35 + 0.47 g from 0.05 mol glucose (1 : 5, 30 min) → 3.7 mmol = 74 mmol/mol [Z, matches
  74]. M 220.4 for C₁₄H₂₄N₂ [Z].
- Preparative isolate from glucose + t-BuNH₂: ~1 mg white residue (product 1), ~0.3 g brown syrup (product 2).
- Spectra: diimine mp 151-152 °C; MS m/e 220 (M⁺, 10), 177 (M - C₃H₇, 74), 82 (100); ¹H-NMR δ 1.30-1.80 (20H),
  3.06 (2H), 7.82 (2H). Glyoxal bis-phenylhydrazone mp 170-171 °C.

## 4. Figures [F] (300 dpi render, ± 5 % of full scale)

Fig. 5 (glucose + t-BuNH₂, 1 M each, 95 % EtOH, reflux):

| t, min | C2 imine (peak area) | 3-DG, mM | radical (ESR) | browning (A420 × dil. / 10) |
|---:|---:|---:|---:|---:|
| 10 | 4.5 | — | 0 | — |
| 15 | — | 0.07 | — | — |
| 20 | 6.7 | — | ~0 | ~0 |
| 30 | 4.7 | — | 0.3 | ~0.1 |
| 35-40 | 4.4 | 1.0 | 1.4 | 0.5 |
| 60 | 3.3 | — | 4.9 | 1.35 |
| 80 | 1.9 | — | 7.4 | 2.0 |
| ~88 | — | 3.6 | — | — |
| 100 | — | — | 7.2 | 2.85 |
| 120 | 1.3 | — | 4.6 | 3.55 |

Fig. 6A x-axis is labelled "t-Butylamine (M)" with ticks to 32 at 1 M glucose — 32 M is impossible for neat
t-butylamine (~9.5 M [Z]); read it as the amine : glucose molar ratio (the text: "maximum at 10 fold excess").
Not two-reader keyed; nothing here is a FIT row.

## 5. What the repository takes

Nothing quantitative for the engine. Context for the oxygen / Namiki item in STATUS_AND_PLAN: (i) an amine-bound
glyoxal-level C2 appears within ~10-20 min at ~78 °C in ethanol and peaks before 3-DG, the radical and colour [F/Q];
(ii) it requires the amine (alkali alone: none) and grows with amine excess to ~10× [Q/F]; (iii) ≥ 1.2-7.4 % of the
glucose can reach it [M]; (iv) pre-Amadori origin from the glycosylamine is proposed, not demonstrated [P].

## 6. Supplementary material

**None referenced.** No supplementary, microfilm or deposited data; the authors defer pathway work to a later paper
("will be published in relation to the formation of the radical species and the browning"). Verdict: nothing held
elsewhere for this paper.

## 7. Flags

1. **Ethanol, alkylamines, no water, no pH** — the main transfer barrier.
2. **Assay reads glyoxal**: the diimine is detected after on-plate hydrolysis; free glyoxal or other glyoxal
   derivatives would be counted with it.
3. **Arbitrary units** for the time course; the yields are isolated (lower bounds) and only for cyclohexylamine.
4. **Fig. 6A axis** mislabelled (see sec. 4).
5. **The oxidation step** implied by Route B is not addressed (no O₂ control).
6. **Registry:** glyoxal dialkylimines are not compound keys (not needed); glyoxal (GO) is.
