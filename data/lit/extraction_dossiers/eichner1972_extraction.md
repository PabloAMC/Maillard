# Eichner & Karel 1972 — EXTRACTION (glucose (or glucose + fructose) + glycine syrups in glycerol-water, a_w 0.23-0.86 set by salt solutions at 37 °C or by glycerol replacing water; ± cellulose, methylcellulose, gum arabic, sorbitol, polyvinyl alcohol; 70 h at 37 / 50 °C or 142 h at 37 °C: ONE end-point browning value (E420 per g glucose) per sample. No Amadori compound, no intermediate, no time course)
DOI: 10.1021/jf60180a025 (not printed; Crossref returns "Influence of water content and water activity on the sugar-amino browning reaction in model systems under various conditions", Eichner, Karel, J. Agric. Food Chem. 20(2) 218-223, March 1972; checked 2026-11-01)

**Read 2026-11-01 from `data/articles/eichner1972.pdf`** (6 pp.; an Acrobat Capture scan). **Page 1 begins with the
reference list of the PRECEDING paper and page 6 continues into the FOLLOWING paper (Donovan et al., ovomucin)**; the
Eichner & Karel article is pp. 218-223 of the journal and is complete. All pages read from the images; **Figs. 1-3 read
from 200 dpi renders** (scratchpad `ing6/eich-2.png`, `eich-3.png`; read error ±0.003 E420 below 0.05, ±0.005
above); Figs. 4-10 read at screen resolution (±0.01). Tables I-III are printed.

**Legend.** **[V]** verbatim · **[M]** printed number · **[F]** read off a figure · **[Q]** qualitative claim · **[C]**
cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Title | "The Influence of Water Content and Water Activity on the Sugar-Amino Browning Reaction in Model Systems Under Various Conditions" |
| Authors | Karl Eichner (present address: Institut für Lebensmitteltechnologie und Verpackung, Munich), Marcus Karel* — Department of Nutrition and Food Science, MIT, Cambridge MA |
| Venue | J. Agr. Food Chem. 1972, **20**(2), 218-223; received 21 July 1971; MIT contribution 1811; US PHS FD-0050 and the German Ministry of Nutrition |
| Lab | MIT (Karel). `lab_ids.yml` has `munster_eichner` for Eichner's later Münster work; by rule 3 a visiting / junior author does not carry a lab, so this is the **MIT Karel group, not in `lab_ids.yml`** (owner decision; Warmbier 1976 / Labuza at Minnesota is a Karel alumnus, separate by the same rule) |

## 1. What this paper answers for the model

The question (atlas sec. 1.2: "Intermediate a_w (0.5 - 0.9) is almost unstudied"): what rates or intermediates vs a_w?

| what the model needs | what Eichner & Karel give | tag |
|---|---|---|
| a rate vs a_w | **no rate**: one E420 after 70 h (or 142 h), called "browning rate" by the authors; no time course, no check of linearity in time | [M] / [!] |
| Amadori / intermediates vs a_w | **none measured** | — |
| browning vs a_w at 37 °C (glucose 10 g + glycine 27.8 mmol, no glycerol; Fig. 2) | E420 **0.083 / 0.191 / 0.172 / 0.072 / 0.014 / 0.005 at a_w 0.32 / 0.41 / 0.52 / 0.62 / 0.75 / 0.84**: maximum at 0.41-0.52; **×38 from 0.84 to 0.41** (×30-50 on the read error), ×14 from 0.84 to 0.62 | [F] / [Z] |
| water content at each a_w (Table III, no glycerol) | **0.09 / 0.15 / 0.20 / 0.31 / 0.50 / 0.86 ml H2O per g glucose** at a_w 0.32-0.84 | [M] |
| is it a_w or water content? | **water content**: at a_w 0.41-0.84 and glycerol 0-0.77 ml/g, all points fall on one curve of E420 vs ml H2O / g glucose (Fig. 3); a_w 0.32 (glassy) falls below it | [F] / [Q] |
| glycerol replacing water at constant reactant concentration (50 °C; Fig. 1) | E420 **0.355 / 0.162 / 0.082 / 0.040 / 0.020** at solvent water 20 / 40 / 60 / 80 / 100 % (a_w ≈0.40 / 0.58 / 0.71 / 0.81 / 0.86, Table II): **×18**, about ×2 per 20 % water removed | [F] / [M] |
| dilution by glycerol at constant water (50 °C; Fig. 5) | at 0.25 ml H2O/g: E420 ≈ 0.49 / 0.32 / 0.165 / 0.08 at 0 / 1 / 2 / 3 ml glycerol/g (×6); at ≥ 1.0 ml H2O/g all ≈ 0.02 | [F] |
| temperature | 37 and 50 °C in different runs (no paired design) | [M] |
| mobility / glass | a_w 0.32: lower than 0.41 "due possibly to the decreased mobility ... very high viscosity" [Q]; glycerol 0.17 ml/g raises it 0.083 → 0.103 (plasticising) [F]; high-viscosity PVA halves browning vs low-viscosity PVA [Q] | [F] / [Q] |

**Answers.**

1. **For the intermediate-a_w regime the paper gives one aggregate response (browning) at 37-50 °C, steeply rising as
   a_w falls from 0.84 to 0.41-0.52 (×14 by 0.62, ×38 by 0.41), with a maximum at 0.41-0.52 for a glucose-glycine
   syrup** [F]. No step-specific rate, no Amadori compound, no dicarbonyl. Confidence that no intermediate is in the
   paper: ~0.99.
2. **The controlling variable is water CONTENT (reactant concentration plus water as a product), not a_w as such**,
   above the glassy range [F/Q]. For a mass-action model this is the useful statement: at fixed temperature, browning
   at 70 h scales roughly as a single function of ml water per g glucose (×~40 from 0.86 to 0.15 ml/g, Fig. 3) [F].
3. **Water as a product inhibitor** (Fig. 1: ×18 at constant reactant concentration when glycerol replaces water) is
   the authors' reading [P]; a solvent effect of glycerol (tautomers, dielectric) is not excluded [Z].

## 2. Conditions

- **Solutions** [M]: (1a-1e) 10 g glucose (55.5 mmol) in 10 ml water or glycerol-water with 20 / 40 / 60 / 80 % glycerol
  (w/w), + 2.085 g glycine (27.8 mmol; partly undissolved in 1e); (2a-2e) the same with 1.043 g glycine (13.9 mmol);
  (3a-3e) 5 g glucose + 5 g fructose. Glucose : glycine 2 : 1 (set 1) or 4 : 1 (set 2) [Z].
- **Runs** [M]: Run 1, set 2, 50 °C, 70 h, closed serum-cap bottles. Run 2, set 1 (1a-1d), 37 °C, 70 h over salt
  solutions (a_w 0.23-0.84, Table I) in equilibrium, water taken up / lost by weight. Run 3, set 2 + extra glycerol,
  50 °C. Run 4, set 2a + 1/3 (w/w) Avicel, Methocel or gum arabic, 37 °C, 70 h, a_w 0.41-0.75. Run 5, sets 2a / 2c /
  2e + 2/3 (w/w) additives incl. glycerol and sorbitol, 50 °C, 70 h, closed. Run 6, set 3 ± 1/3 PVA (Elvanol 50-42 or
  51-05), 37 °C, **142 h**, a_w 0.23-0.62.
- **pH**: not stated, not controlled, no buffer [!]. O2: air in serum-cap bottles or desiccators, not controlled.
- **Analysis** [M]: diluted to E420 0.1-1.0, reported per 1 g glucose in 100 ml; unheated mixtures as blanks (Run 4).
  a_w of the glycerol solutions by weight change over salts (graphical interpolation). Replicates: not stated.

## 3. Quantitative results

**Table I** [M]: salt solutions a_w at 37 °C: CH3COOK 0.23, MgCl2 0.32, K2CO3 0.41, Mg(NO3)2 0.52, NaNO2 0.62, NaCl
0.75, K2CrO4 0.84. **Table II** [M]: a_w of 2a-2e at 37 °C: 0.86 / 0.81 / 0.71 / 0.58 / ≅0.40.

**Table III** [M] (Run 2: ml H2O / g glucose at glycerol 0 / 0.17 / 0.35 / 0.55 / 0.77 ml/g):

| a_w | 0 | 0.17 | 0.35 | 0.55 | 0.77 |
|---:|---:|---:|---:|---:|---:|
| 0.32 | 0.09 | 0.13 | 0.17 | 0.20 | 0.22 |
| 0.41 | 0.15 | 0.19 | 0.24 | 0.30 | 0.37 |
| 0.52 | 0.20 | 0.25 | 0.31 | 0.35 | 0.41 |
| 0.62 | 0.31 | 0.40 | 0.47 | 0.55 | 0.63 |
| 0.75 | 0.50 | 0.62 | 0.78 | 0.93 | 1.10 |
| 0.84 | 0.86 | 1.11 | 1.37 | 1.73 | 2.07 |

(Read from a 200 dpi render of Table III; the OCR layer misreads the a_w 0.32 / 1e glycerol cell as "0.17", the image
prints 0.77.)

**Fig. 2** [F] (Run 2, 37 °C, E420 vs glycerol 0 / 0.17 / 0.35 / 0.56 ml/g): a_w 0.32: 0.083 / 0.103 / 0.072 / 0.059;
0.41: 0.191 / 0.119 / 0.082 / 0.056; 0.52: 0.172 / 0.106 / 0.058 / 0.020; 0.62: 0.072 / 0.039 / 0.022 / 0.012; 0.75:
0.014 / 0.009 / 0.006 / 0.005; 0.84: 0.005 / 0.003 / 0.002 / 0.002.

**Fig. 1** [F] (Run 1, 50 °C): see sec. 1. **Fig. 3** [F]: one curve, E420 0.19 at 0.15 ml/g → 0.05 at ~0.33 → 0.01
at ~0.6 → ~0.003 at ≥ 0.9 ml/g.

**Run 4 / 5 (Figs. 6-7)** [F]: Methocel raises browning (e.g. a_w 0.41: 0.17 vs 0.10 original), gum arabic lowers it
at low a_w (0.04 at 0.41); Avicel ≈ original; at 20 % solvent water (50 °C) original 0.53, Methocel 0.60, Avicel 0.42,
gum arabic 0.28, glycerol 0.27, sorbitol 0.14.

**Run 6 (Figs. 8-10)** [F] (glucose + fructose, 142 h, 37 °C): without PVA the maximum is at a_w 0.32-0.41 (E420 ≈
0.05-0.06); glycerol first raises then lowers; high-viscosity PVA gives ~50 % of the low-viscosity PVA browning [Q].

## 4. Mechanism evidence and its kind

- **Within-study factors** (kind 6): a_w at fixed composition; water content at fixed a_w; glycerol dilution at fixed
  water; viscosity (PVA grade). Single end-point, no fit, no intermediate.
- Mechanism claims are [P]: water as a condensation product inhibits by mass action; low mobility in viscous / glassy
  syrups lowers the rate at a_w ≤ 0.32.
- Lab: MIT (Karel), with Eichner visiting. One paper.

## 5. Bearing on open questions

- **The intermediate-a_w gap (atlas sec. 1.2; evidence map 1.6, "Deoxyosones on an a_w ladder: NO study"):
  LEAVES UNCHANGED for intermediates** (none measured). For the aggregate response it adds a second
  concentration-framed data set beside Labuza / Warmbier and Roos: browning is set by water content (mass action +
  product water), which a mass-action engine can express if it carries water as a species and solute concentrations
  rise as water falls [Z].
- **F05 / F31 (a_w tags):** SUPPORTS "a_w maximum, system-dependent": 0.41-0.52 for a protein-free glucose-glycine
  syrup at 37 °C; Warmbier 1976 (glycerol, casein) puts lysine loss's maximum at 0.45-0.55
  (`warmbier1976_extraction.md`) — consistent. Glycerol plasticises below the maximum and dilutes above it, the same
  two-sided effect Warmbier reports.
- **E5 / E7 / CF-011 / CF-031:** not touched.

## 6. Disagreements found

1. **None with another dossier on the numbers.** The paper's statement that the maximum is NOT fixed by a_w is the
   point: any engine law that keys browning to a_w alone (rather than water content and mobility) is contradicted
   here [Q].
2. **[!] Run 2 a_w for solution 1e** (undissolved glycine) is omitted from Fig. 2; Table III still lists 1e. Fig. 2's
   abscissa stops at 0.56 ml/g, so the 0.77 ml/g column is not plotted.
3. **[!] "Browning rate"** is the authors' word for a 70 h end-point; it is a rate only if E420 is linear in time over
   70 h, which is not shown. Treat every value as an extent.

## 7. Ledger-ready measurements

1. **Fig. 2 at glycerol 0** (two readers): E420 at six a_w, 37 °C, with the water content of Table III (printed) — the
   a_w / water-content pair for a glucose-glycine syrup.
2. **Fig. 1** (two readers): five glycerol / water ratios at constant reactant concentration, 50 °C (a_w from Table II)
   — the water-as-product test.
3. Fig. 3 (two readers): the collapse onto water content, as a check on any a_w law.
