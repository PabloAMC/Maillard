# Wolf, Thompson & Reineccius 1977 — EXTRACTION
DOI: 10.1111/j.1745-4549.1977.tb00330.x (not printed in the PDF; Crossref 2026-11-01: "KINETICS OF AVAILABLE LYSINE LOSS DURING THERMAL PROCESSING OF A SOY PROTEIN ISOLATE: REVERSION OF NO-LOSS TO RAPID-LOSS PHASE", Wolf, Thompson, Reineccius, J. Food Process. Preserv. 1(4), 271-278, 1977-10. **Matches the brief's DOI.**)
### Soy protein isolate (Promine D) + microcrystalline cellulose +- 4 % glucose (dry basis), slurried in 0.1 M pH 7.0 sodium phosphate, freeze-dried, equilibrated to a_w 0.65 (saturated CuCl2), sealed in retort pouches, 125 °C for 2 min (lag-corrected), FDNB available lysine. Two experiments: heat the protein without sugar first, then add glucose; and shear / pressure between two processings. **No rate constant, no Ea, one temperature, one a_w.** The kinetic constants the brief asks for are in the companion papers (Wolf et al. 1977a, J. Food Sci. 42:1540; 1978, J. Food Sci. 43:1486), not on disk.

**Read 2026-11-01 from `data/articles/wolf1977.pdf`** (8 pp., scanned journal pages with an OCR layer). Prose and the
two tables from `pdftotext -layout`; **Fig. 1 (schematic retention curve) and Table 1 checked against 150 dpi renders
of pp. 2 and 5** (scratchpad `ing5/wolf-2.png`, `wolf-5.png`). Fig. 1 is a drawn "typical" curve without data points;
read error +-0.3 mg/g, +-0.5 min.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** University of Minnesota, St Paul: Agricultural Engineering (Wolf, Thompson) and Food Science
(Reineccius). **Existing lab id `minnesota_reineccius`** (not `minnesota_labuza`); not new.

## 1. What this paper answers for the model

| model need | what Wolf 1977 gives | tag |
|---|---|---|
| available-lysine loss rate constant(s) | **none printed.** The intro says loss is first-order, then stops after **33-55 %** loss (read from the p. 271 render; the OCR layer says "3345%") [M/C], and cites Wolf 1977a / b for the kinetics | [C] |
| Ea | **none** (one temperature, 125 °C; Fig. 1 schematic at 105 °C) | — |
| a_w / moisture | a_w 0.65 set by saturated CuCl2 in a partial vacuum (not instrument-verified); moisture content not printed; sealed retort pouches, so moisture is held during heating | [M] |
| pH | slurry 0.1 M sodium phosphate pH 7.0 before freeze-drying; Fig. 3 re-adjusts to pH 7.0 before re-processing; in-situ pH not measured | [M] |
| glucose : lysine | 4 % glucose dry basis; lysine 9.6-10.2 mg/g sample -> glucose : lysine ~3.4 mol/mol [Z] | [M] / [Z] |
| loss in one 125 °C / 2 min pass | **9.6 -> 5.0 mg/g (-48 %)** with glucose; **10.2 -> 4.8-5.3** in the shear / pressure set | [M] |
| second pass (re-slurried, re-dried, a_w 0.65 again) | **no further loss**: 5.3 -> 4.8 (shear), 5.3 -> 5.3 (pressure), 4.8 -> 5.1 (neither); F = 1.68 / 0.08 / 3.63, n.s. | [M] |
| protein heated WITHOUT sugar, then glucose added | 4.8 +- 0.5 mg/g, the same as glucose present from the start (5.0 +- 0.7): heat alone gives no protection | [M] |

**Answers.**
1. **Rate constants and Ea: not in this paper** (confidence ~97 %). What it contains is a two-point loss per pass. A
   lower bound on the first-pass pseudo-first-order constant at 125 °C, a_w 0.65 [Z]: ln(9.6 / 5.0) / 2 min =
   **>= 0.33 min⁻¹** (a bound: the curve may already have entered its no-loss phase inside the 2 min, and the lag
   correction is not printed). Fig. 1 (schematic, 105 °C) [F]: ~15 -> ~6.7 mg/g in ~11 min, i.e. ~0.07 min⁻¹ [Z], a
   drawn curve, not data.
2. **Against the engine's glycation step (F05)** [Z, transfer-laden]. Frozen B20 k_glyc = 10^-4.707 L/(mmol min) at
   100 °C with Ea 96.8 kJ/mol gives **1.39 × 10⁻⁴ L/(mmol min) at 125 °C**. Matching 0.33 min⁻¹ by plain mass action
   needs ~2.3 M glucose in the water phase, i.e. 0.222 mmol glucose/g dry solids in <= ~9.5 g water/100 g solids. The
   moisture at a_w 0.65 is not printed (soy-cellulose sorption at a_w 0.65 is typically ~8-12 g/100 g; NOT from this
   paper). So **the first pass is consistent with the engine's constant within ~x2**, given the water-phase
   concentration convention. **The second pass is not**: with glucose still in 3x molar excess the engine would remove
   another ~40-50 % of the lysine on re-processing; Wolf finds none (Table 2). The engine has no "no-loss phase"
   (an inaccessible / protected lysine pool). Ea cannot be tested.
3. **CF-072 (does the bound Amadori compound return lysine?)** — the schematic Fig. 1 shows lysine RISING in a
   "transition phase" after the minimum (~6.7 -> ~9.3 mg/g) [F, schematic]. Read literally it is a partial lysine return,
   the Brands-type structure the engine writes. Weight: **low (~0.15)**: it is a drawn "typical" curve; the authors
   interpret it as a change of kinetics, not a return; and FDNB "available" lysine is known to count part of the
   early Maillard lysine as available (my background knowledge, not from this paper) [Z].

## 2. Conditions

- **System** [M]: "20% Protein (Promine D, Central Soya)", glucose, microcrystalline cellulose (Avicel PH-101); glucose
  added before freeze-drying (Fig. 3) or after a first sugar-free processing (Fig. 2).
- **Buffer / pH** [M]: dry-blended, slurried with 0.1 M pH 7.0 sodium phosphate (v/w 1.5 : 1), frozen at -30 °C,
  freeze-dried, crushed.
- **Water** [M]: equilibrated to a_w 0.65 in desiccators over saturated cupric chloride (partial vacuum).
- **Heating** [M]: sealed retortable foil pouches, miniature retorts, 125 °C for 2 min; heating / cooling lags "offset"
  by a Gondo-type correction (the corrected equivalent time is not printed).
- **Shear / pressure** [M]: Omni-mixer 10,000 rpm 2 min (> 1000 s⁻¹, ice bath); homogeniser > 500 psi (3.45 MPa); then
  pH 7.0, freeze-dry, a_w 0.65, re-process 125 °C / 2 min.
- **Analysis** [M]: Carpenter FDNB available lysine (Booth modification); glucose by hexokinase kit (**no glucose value
  is printed** [!]); 4 replicates; ANOVA with "week of analysis" as covariate.

## 3. Quantitative results

**Table 1** (available lysine, mg/g sample, mean +- SD, n = 4) [M]: 1 unprocessed, no glucose 9.6 +- 1.2; 2 processed
once, no glucose 7.5 +- 1.1; 3 twice, no glucose 9.1 +- 1.2; 4 twice, glucose added before the second 4.8 +- 0.5;
5 once with 4 % glucose 5.0 +- 0.7 (from Wolf 1977a).

**Table 2** (available lysine remaining, mg/g sample) [M]: control (no process) 10.2 +- 1.0; shear 5.3 +- 0.7 -> 4.8 +-
0.4 (F 1.68); pressure 5.3 +- 0.4 -> 5.3 +- 0.3 (F 0.08); no shear or pressure 4.8 +- 0.4 -> 5.1 +- 0.2 (F 3.63); pooled
F 0.95; all n.s. at 0.05.

**Fig. 1** (schematic, 105 °C) [F]: ~15 mg/g at t = 0, minimum ~6.7 at ~11 min ("1st order rapid loss phase"), rise to
~9.3 by ~17 min ("transition phase"), flat to 35 min ("no loss phase").

## 4. Mechanism evidence and its kind

Within-study factors (sugar-free pre-heating; shear; pressure; re-processing) on one marker. Hypotheses eliminated
[Q]: sugar exhaustion and local sugar depletion (by Wolf 1977a / b, [C]), heat-induced protection without sugar
(Table 1), reversion by shear or pressure alone (Table 2). Mechanism of the no-loss phase: unidentified, stated as
such. Laboratory: Minnesota (Reineccius).

## 5. Bearing on open questions

- **F05 unknown "pea / soy isolate glycation KINETICS at cooking temperature"** — PARTLY: a soy-isolate loss at 125 °C
  (>= 0.33 min⁻¹ first pass, then a plateau), not a fitted constant. The constants exist in Wolf 1977a / 1978 (lead
  below).
- **F05 unknown "high T (100-140 C) at intermediate a_w with a measured a_w"** — PARTLY filled: 125 °C at a_w 0.65 in
  SEALED pouches (moisture held), a_w set by a saturated salt but not instrument-verified; one temperature.
- **Atlas sec. 1.4 "plant-protein-bound amines"** — PARTLY filled (soy, FDNB lysine, two points per pass).
- **F05 engine mismatch (d), no a_w term on glycation** — LEAVES UNCHANGED (one a_w).
- **CF-072** — LEAVES UNCHANGED in substance; the schematic's lysine rise is a weak qualitative hint toward a return
  step (sec. 1, answer 3).
- **Engine structure** — CONTRADICTS a pure mass-action glycation with ample sugar (no plateau in the engine) at
  a_w 0.65, 125 °C; see sec. 6.

## 6. Disagreements found

1. **Engine vs Wolf: lysine loss stops at ~50 % with sugar in excess (NEW).** Quantity: available-lysine loss in a
   second 125 °C / 2 min pass at a_w 0.65. Wolf: 0 (-0.5 to +0.3 mg/g, n.s.); engine (k_glyc as frozen, glucose 3.4x
   molar excess, water-phase concentration): another ~40-50 % [Z]. Candidate causes: a protected / inaccessible
   lysine pool after the first glycation (authors' hypothesis space: conformation, but only with glucose + heat);
   loss of glucose to non-lysine sinks during the first pass (glucose not printed here, asserted adequate in Wolf
   1977a [C]); FDNB counting Amadori lysine as available (would flatten apparent loss); the re-slurry / re-dry
   cycle (pH reset to 7.0) changing the matrix.
2. **Unexplained sugar-free "loss"** [!]: Table 1 treatment 2 (processed once, no glucose) 7.5 vs 9.6 mg/g untreated
   (-22 %, ~1.3 pooled SD), but treatment 3 (processed twice, no glucose) 9.1. Read as noise; the authors do not
   discuss it.

## 7. Ledger-ready measurements

- Table 2 and Table 1 rows 1, 4, 5: available lysine before and after one and two 125 °C / 2 min passes at a_w 0.65,
  pH 7.0 (pre-dry), 4 % glucose — a soy-isolate end-point pair showing the plateau. Key as end points (two readers),
  with the lag correction unstated.
- Do NOT key Fig. 1 (schematic).
- **Download lead** (Crossref-checked 2026-11-01): Wolf, Thompson & Reineccius 1977a, "Initial losses of available
  lysine in model systems", J. Food Sci. 42:1540-1544, 10.1111/j.1365-2621.1977.tb08420.x — the first-order constants
  and conditions; and 1978, "Comparisons between model predictions and measured values for available lysine losses in
  a model food system", J. Food Sci. 43:1486-1490, 10.1111/j.1365-2621.1978.tb02525.x.
