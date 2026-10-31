# Warmbier, Schnickels & Labuza 1976 — EXTRACTION (solid intermediate-moisture model food: casein 30 g + glucose 10 g + GLYCEROL 20 g + oil + cellulose, glucose : available lysine 2.9 mol/mol; 25 / 35 / 45 °C, a_w 0.14 - 0.87, sealed cans, up to ~100 d; INITIAL first-order loss of glucose and FDNB-available lysine as RATES (mmol / 100 g solid / day) vs a_w, maximum at a_w ~0.45 - 0.55; brown pigment zero order after an induction time, Ea 32.9 kcal/mol; ~1.3 mol glucose lost per mol lysine; up to ~70 % of the lysine gone before visible browning)

### A PROTEIN-BOUND-LYSINE DATA SET WITH A HUMECTANT. The a_w maximum for lysine loss sits at 0.45 - 0.55 — LOWER than the 0.65 - 0.75 of glycerol-free powders — because liquid glycerol keeps reactants mobile below the water monolayer. Rates are printed only as figure points (initial rates), not as rate constants; the rate constants follow from the formulation [Z].

**Source on disk:** `data/articles/warmbier1976.pdf` (4 pp., J. Food Sci. 41:528-531; Acrobat Capture OCR of a scan).
Prose and Tables 1-2 from the `pdftotext -layout` layer, checked against **200 dpi renders**, from which **Figs. 1-5
were read** (scratchpad `pdf3/war-2.png`, `war-3.png`).

**Legend.** **[V]** verbatim (short quote) · **[M]** a number printed in this paper · **[F]** read from a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed / hypothesised, not
measured · **[Z]** derived by me, never printed.

## 1. What this paper answers for the model

The question (`docs/guides/MECHANISM_REVIEWS_2026-10-31.md` sec. 3 Tier 1; `docs/history/mechanism_reviews_2026-10-31/
matrix.md` gaps (a) bound lysine and (b) low a_w): kinetics of available-lysine loss in casein + glucose vs water
activity — rate constants, order, temperature, a_w dependence, glucose : lysine stoichiometry, browning. The engine's
a_w multiplier (`src/kinetic_core/trunk_conditions.py`, `AW_STEPS = ("k_schiff",)`) comes from Pereyra Gonzales 2010
(skim milk powder, no humectant), held flat below a_w 0.33.

| what the model needs | what Warmbier 1976 gives | tag |
|---|---|---|
| order | "The initial loss rate of both glucose and available lysine, however, follows first order kinetics" (abstract); Fig. 3 lysine curves decelerate after ~1 week (initial only) | [V] / [F] |
| lysine loss rates vs a_w | **initial rates, mmol lysine / 100 g solid / day**: 45 °C — ~1.5 (a_w 0.29), ~2.1 (0.47), **~2.25 (0.57, max)**, ~1.05 (0.78), ~0.85-0.97 (0.85-0.86); 35 °C — ~0.18-0.30 (0.33-0.34), ~0.32-0.46 (0.41-0.42), **~0.9-1.0 (0.50-0.52, max)**, ~0.44-0.48 (0.61-0.65), ~0.65-0.70 (0.73), ~0.40 (0.79), ~0.22 (0.87). DM and HM water addition overlap | [F] |
| lysine first-order k | **k ≈ 0.12 d⁻¹ at 45 °C, a_w ~0.57; ≈ 0.05 d⁻¹ at 45 °C, a_w ~0.85; ≈ 0.05 d⁻¹ at 35 °C, a_w ~0.5; ≈ 0.012 d⁻¹ at 35 °C, a_w 0.87** (rate ÷ 19.1 mmol available lysine / 100 g); Fig. 3 directly: ~0.15-0.17 d⁻¹ (45 °C, a_w 0.57), ~0.07 d⁻¹ (a_w 0.85) | [Z] |
| glucose first-order k | 45 °C peak ~3.7 mmol / 100 g / d at a_w ~0.47 → **0.067 d⁻¹**; Fig. 3 (a_w 0.56) gives 0.068 d⁻¹ — the basis check that "100 g solid" = 100 g of the water-free formulation; 35 °C peak ~1.03 → 0.019 d⁻¹ | [F] / [Z] |
| temperature | pigment: Ea **32.9 kcal/mol** (137.7 kJ/mol [Z]), Q10 6.1 (25-35 °C) and 5.4 (35-45 °C), 33× from 25 to 45 °C [M]. No Ea for lysine or glucose; from the a_w-maximum rates 35 → 45 °C: lysine ×~2.25 → **Ea ~66 kJ/mol**, glucose ×~3.6 → **~104 kJ/mol** [Z, two points, maxima at different a_w] | [M] / [Z] |
| a_w dependence | bell-shaped; maximum for pigment at a_w 0.45-0.55, "maximum loss rate of these reactants occurs near the same a_w" [Q]; at the maximum, moisture 8.9 (35 °C) / 10 (45 °C) g H₂O / 100 g solids vs BET monolayer 8.1 [M]. Ratio max / a_w ~0.86: lysine ×2.3-2.6 (45 °C), ×4.5 (35 °C) [Z] | [M] / [F] / [Z] |
| glucose : lysine stoichiometry | initial charge 2.9 mol glucose / mol available lysine [M]; "an average of 1.3 moles of glucose initially react for one mole of lysine" at 35 and 45 °C [V]; Table 2 end-of-induction ratios scatter 0.7-4.3 [Z] | [M] / [V] / [Z] |
| browning | zero order after an induction time of **100 d (25 °C), 20 d (35 °C), 7 d (45 °C)**; peak pigment rates ~0.012 (35 °C, a_w ~0.4) and ~0.051 (45 °C, a_w ~0.55) ΔA420 / g solid / day; 10-33 % glucose and 12-71 % lysine lost during the induction (Table 2) | [M] / [F] |

**Answers.**

1. **Rate constants: derivable, not printed** (confidence in the derivation ~80 %: it rests on reading figure points
   to ± 10 % and on the "100 g solid" basis, which the glucose cross-check supports). Lysine k ≈ 0.12 d⁻¹ (45 °C) and
   ≈ 0.05 d⁻¹ (35 °C) at the a_w maximum — pseudo-first-order in a 2.9-fold glucose excess, so they fold in the
   glucose level and should NOT be read as per-amine second-order constants without a concentration basis (the
   aqueous phase volume is not given; glycerol is part of the solvent).
2. **The a_w law is matrix-specific** (P ~0.8): with liquid glycerol the maximum is at a_w 0.45-0.55 and the high-a_w
   decline is ×2.3-4.5 by a_w 0.85-0.87; the authors attribute the downward shift to glycerol's "water-like
   properties" [Q, P]. The engine's single skim-milk curve cannot serve a glycerol-humectant product (e.g. the
   matrix review's glycerol-adjusted a_w design) and a glycerol-free powder at once.
3. **Stoichiometry > 1**: 1.3 mol glucose per mol lysine lost initially — glucose has sinks besides the lysine
   adduct (or multiple glycation) [V/Z]; the engine's 1 : 1 Schiff step plus sugar degradation is consistent in
   sign only.
4. **Nutritional loss precedes colour**: up to 69-71 % of the available lysine is gone at 45 °C before zero-order
   browning starts (Table 2) — an early-stage (Amadori) sink, not melanoidin.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of Glycerol on Nonenzymatic Browning in a Solid Intermediate Moisture Model Food System" |
| Authors | H. C. Warmbier, R. A. Schnickels, T. P. Labuza — Dept. of Food Science & Nutrition, University of Minnesota, St. Paul |
| Venue | *J. Food Sci.* 1976, **41**(3), 528-531; MS received 7/11/75, revised 10/3/75, accepted 10/6/75 |
| DOI | 10.1111/j.1365-2621.1976.tb00663.x (not printed; Crossref returns this title, these three authors, J. Food Sci. 41(3) 528-531, 1976; checked 2026-10-31). **Matches the brief.** |
| Funding | Univ. Minnesota Agric. Expt. Stn. Project 18-72; NASA contract NAS 9-12560 (Johnson Space Center) |
| Paper type | Storage kinetics; FDNB-available lysine (Booth 1971), glucose oxidase kit, trypsin-digest A420 pigment |
| Related | Warmbier, Schnickels & Labuza 1976, J. Food Sci. 41:981 (glucose : lysine ratio; listed in `matrix.md`, not on disk); Warmbier 1975 PhD thesis (Univ. Minnesota) — holds the 25 °C a_w series; Eichner & Karel 1972 JAFC 20:218 (glycine-glucose-glycerol, maximum a_w 0.41 at 37 °C) [C] |

## 2. Methods as they matter to a model

- **Formulation** (Table 1) [M]: K-sorbate 0.3 g, glucose 10.0, glycerol 20.0, casein 30.0, Apiezon B oil 20.0,
  microcrystalline cellulose 20.0, water variable. Water-free total 100.3 g [Z]. Glucose 55.3 mmol / 100 g; available
  lysine 55.3 / 2.9 = **19.1 mmol / 100 g** = 0.64 mmol / g casein [Z].
- **Water** [M]: Direct Mix (liquid water, 1 d over saturated salt) vs Humidified Mix (3-6 d over saturated salt);
  moisture by methanol-extraction GC, a_w by vapour-pressure manometry. a_w 0.14-0.87.
- **Storage** [M]: 202 × 214 epoxy-lined cans, sealed, ends dipped in glyptol; 25, 35, 45 °C; held at -29 °C if
  analysis delayed. No pH, no headspace / O₂ control reported.

## 3. Numbers printed

Table 2 (loss during the browning induction time) [M]:

| T, °C | a_w | induction, d | glucose % loss | lysine % loss |
|---:|---:|---:|---:|---:|
| 25 | 0.86 | 100 | 10 | 60 |
| 25 | 0.72 | 100 | 21 | 65 |
| 25 | 0.67 | 100 | 26 | 55 |
| 35 | 0.87 | 20 | 13 | 20 |
| 35 | 0.51 | 20 | 33 | 60 |
| 35 | 0.32 | 20 | 18 | 12 |
| 45 | 0.84 | 7 | 29 | 37 |
| 45 | 0.39 | 7 | 33 | 69 |
| 45 | 0.14 | 7 | 17 | 71 |

Text [M]: pigment Ea 32.9 kcal/mol, Q10 6.1 / 5.4, ×33 (25 → 45 °C); maximum pigment rate a_w 0.45-0.55; moisture at
maximum 8.9 / 10 g H₂O / 100 g solids; BET monolayer 8.1; 1.3 mol glucose / mol lysine.

## 4. Figures [F] (200 dpi render; ± 10 % of the value)

- Fig. 1: A420 / g solids vs days; 45 °C (a_w 0.73) reaches ~0.8 by ~17 d; 35 °C series (a_w 0.51, 0.73, 0.87) reach
  ~0.43 / ~0.38 / ~0.2 at ~50 d; 25 °C (a_w 0.74) ~0.02 at 53 d.
- Fig. 2: ΔA420 / g solid / day vs a_w; 35 °C peak ~0.0118 at a_w ~0.4, ~0.0046 at 0.9; 45 °C peak ~0.051 at a_w
  ~0.55, ~0.016 at 0.15, ~0.026 at 0.86. **Caption says "25 and 35 °C" but the panels are labelled 35 and 45 °C.**
- Fig. 3 (45 °C, DM): fraction remaining; glucose (log axis) a_w 0.56: ~0.86 (2 d), ~0.64 (6.5 d), ~0.48 (21 d);
  a_w 0.73: ~0.55 (21 d); a_w 0.88: ~0.75 (21 d). Lysine (axis labelled 0-1.0, near-linear) a_w 0.57: ~0.54 (3.7 d),
  ~0.38 (6.5 d), ~0.32 (15 d); a_w 0.78: ~0.64 (6.5 d), ~0.48 (15 d); a_w 0.85: ~0.62 (6.5 d), ~0.44 (15 d).
- Figs. 4-5: initial rates as tabulated in sec. 1.

Not two-reader keyed. If any number here becomes a FIT row, key Figs. 3-5 by two readers first.

## 5. What the repository takes

Nothing is wired yet. Candidates for the matrix / a_w work (roadmap step 6), not priors: (i) protein-bound (casein)
lysine loss k ≈ 0.05-0.12 d⁻¹ at 35-45 °C at the a_w maximum, ×2-4.5 lower at a_w 0.85-0.87 [Z]; (ii) a
glycerol-humectant a_w curve peaking at 0.45-0.55 as a rival to the skim-milk curve for glycerol products [F/Q];
(iii) 1.3 mol glucose per mol lysine lost [V]; (iv) pigment Ea 137.7 kJ/mol with a T-dependent induction time [M].

## 6. Supplementary material

**No supplementary, microfilm or deposited data.** The paper points to **Warmbier 1975, PhD thesis, University of
Minnesota** for the 25 °C a_w series ("The same pattern was also found at 25 °C") — the 25 °C rates and the raw
time courses the model would need for real rate constants live there, not here. The companion J. Food Sci. 41:981
(glucose : lysine ratio) is a separate paper. Verdict: the thesis holds data the model needs (25 °C; full curves);
low priority unless the a_w law is fitted.

## 7. Flags

1. **Rates, not constants**; figure-read; "100 g solid" basis inferred (glucose cross-check agrees within ~2 %).
2. **"Initial" first order only**: Fig. 3 lysine curves flatten after ~1 week (a non-reactive or blocked lysine
   fraction, or reversibility of the FDNB-reactive adduct) — a single first-order law over-predicts late loss.
3. **FDNB-available lysine** counts the ε-amine; Amadori-bound lysine partly regenerates under FDNB hydrolysis in
   some protocols (not discussed here) [Z].
4. **No pH, no O₂ / headspace data**; casein type not specified.
5. **Fig. 2 caption temperatures wrong** (see sec. 4).
6. **Registry:** casein / protein-bound lysine has no compound key yet (needed only when roadmap step 6 lands).
