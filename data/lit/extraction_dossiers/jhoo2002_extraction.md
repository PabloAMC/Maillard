# Jhoo et al. 2002 — EXTRACTION (2-methyl-3-furanthiol trapped by thiamin as the thioether MAMP)

**Source on disk:** `data/articles/jhoo2002.pdf` (the publisher's PDF, 4 pages, pp. 4055-4058). Read and
checked by eye on 2026-10-09: every number below was read from the page image and cross-checked against
`pdftotext -layout`. Tables 2 and 3 are set as images with drawn structures and have no text layer at
all, so their numbers were read from a 300 dpi render of p. 4057 only; the compound names in those two
tables are not printed, only structures, and the names used below are assigned here from the drawings.
Written for the thiamine route of the sulfur lane (`src/kinetic_core/sulfur.py`, "THE THIAMINE ROUTE"),
to answer one question: does thiamin consume MFT?

| field | value |
|---|---|
| Title | "Characterization of 2-Methyl-4-amino-5-(2-methyl-3-furylthiomethyl)pyrimidine from Thermal Degradation of Thiamin" |
| Authors | Jin-Woo Jhoo, Ming-Chi Lin, Shengmin Sang, Xiaofang Cheng, Nanqun Zhu, Ruth E. Stark, Chi-Tang Ho (Rutgers; College of Staten Island, CUNY) |
| Venue | J. Agric. Food Chem. 2002, 50(14), 4055-4058. Received 4 Dec 2001, revised 26 Apr 2002, accepted 29 Apr 2002, published on the Web 05/31/2002. |
| DOI | 10.1021/jf011591v (printed, p. 4055 foot) |

MAMP = 2-methyl-4-amino-5-(2-methyl-3-furylthiomethyl)pyrimidine, compound **1**; MFT = 2-methyl-3-
furanthiol, compound **7**; the thiazole 4-methyl-5-(2-hydroxyethyl)thiazole is compound **6**.

## 1. Methods (pp. 4055-4056)

Four experiments, at different conditions. They must not be pooled.

| experiment | charge | solvent | heating | analysis |
|---|---|---|---|---|
| Preparative degradation (p. 4055) | thiamin hydrochloride 30 g | 0.1 M potassium phosphate buffer, 300 mL, **pH 6.5** | 1000 mL flask with condenser, **"at 110 °C in an oil bath for 2 h"** | CH₂Cl₂ extraction (3 × 300 mL), silica column, RP-18; isolation and NMR/MS structure proof |
| Model system I (p. 4056) | thiamin monochloride 1 g | **methanol** 15 mL | 80 °C, 2 h, oil bath | CH₂Cl₂ extraction, GC-MS (no internal standard stated) |
| Model system II (p. 4056) | thiamin monochloride 1 g + MFT 1 g | **methanol** 15 mL | 80 °C, 2 h, oil bath | as I |
| Model system III (p. 4056) | thiamin monochloride 1 g | **water** 15 mL | "refluxed in an oil bath at 80 °C for 70 min" | tridecane internal standard (1 mL of 1000 ppm), CH₂Cl₂ 15 mL, GC-MS |
| Model system IV (p. 4056) | thiamin monochloride 1 g + cysteine 1 g | **water** 15 mL | as III | as III |

Also: MAMP synthesis for the structure proof, thiamin monochloride 1 g + MFT 1 g in methanol 15 mL,
80 °C, 2 h (p. 4056); thiamin monochloride itself prepared from thiamin HCl with triethylamine in
methanol.

Notes on the conditions, as printed:

- **pH is stated only for the preparative run (6.5).** Not for I-IV.
- The 110 °C is the **oil-bath** temperature for an aqueous solution under a condenser; the solution
  temperature is not printed (it cannot have been far above its boiling point — inference, not stated).
- **III and IV contradict themselves:** the methods say water (p. 4056); the discussion describes III
  and IV as "thiamin monochloride in methanol" with and without cysteine (p. 4058). The methods paragraph is the more specific and is taken here; the conflict is unresolved.
- Concentrations, derived here (thiamin monochloride C₁₂H₁₇ClN₄OS 300.81 g/mol; thiamin HCl 337.26;
  MFT 114.16; cysteine 121.15):
  preparative 30 g / 337.26 = 89.0 mmol in 0.300 L = **0.297 M**;
  I-IV thiamin 1 g / 300.81 = 3.32 mmol in 0.015 L = **0.222 M**;
  II MFT 1 g / 114.16 = 8.76 mmol = **0.584 M**; IV cysteine 1 g / 121.15 = 8.25 mmol = **0.550 M**.

## 2. Findings that matter

### Preparative run (p. 4056): MAMP is a real product of thiamin alone in pH 6.5 buffer

CH₂Cl₂ extract 700 mg → MAMP 20 mg isolated, C₁₁H₁₃ON₃S, HRFAB-MS [M+H]⁺ m/z 236.0865 (calc.
236.0858). Structure from ¹H/¹³C/2D NMR (Table 1, p. 4056; the bridging CH₂ carbon C-8 at δ 34.3 in
MAMP vs 51.0 in thiamin) and confirmed by synthesis. The paper calls it "a major decomposition product".
*Derived here:* 20 mg / 235.30 g/mol = 0.085 mmol from 89.0 mmol thiamin HCl = **0.096 mol % isolated
yield** (a lower bound: isolation losses unstated); 20/700 = 2.9 % of the extract mass. No MFT
quantity is reported for this run.

### Table 2 (p. 4057): concentration, mg/mol of thiamin, methanol, 80 °C, 2 h

| compound (drawn structure; name assigned here) | Model system I (thiamin) | Model system II (thiamin + MFT) |
|---|---|---|
| 2-methyl-3-furanthiol (MFT, **7**) | — (blank) | 211 |
| 4-methyl-5-(2-hydroxyethyl)thiazole (**6**) | 1147 | 163 |
| 4-amino-5-(methoxymethyl)-2-methylpyrimidine | 307 | — (blank) |
| bis(2-methyl-3-furyl) disulfide | — (blank) | 85 |
| MAMP (**1**) | 32 | 211 |

Footnote a: "(1) thiamin monochloride; (2) thiamin monochloride and 2-methyl-3-furanthol" [sic].

### Table 3 (p. 4057): concentration, mg/mol of thiamin, 80 °C, 70 min (water per the methods)

| compound (drawn structure; name assigned here) | Model system III (thiamin) | Model system IV (thiamin + cysteine) |
|---|---|---|
| 2-methyl-3-furanthiol (MFT, **7**) | — (blank) | 11 |
| thiophene-3-thiol (drawn: thiophene, SH at C-3, no methyl) | — (blank) | 134 |
| 4-methyl-5-(2-hydroxyethyl)thiazole (**6**) | 1474 | 12894 |
| MAMP (**1**) | 29 | 12 |

Footnote a: "(1) thiamin monochloride; (2) thiamin monochloride and cysteine." Blank cells are printed
blank; no detection limit is given, so a blank means "not reported", not a measured zero.

### Molar yields — derived here (mg/mol ÷ g/mol = mmol per mol thiamin)

| | I | II | III | IV |
|---|---|---|---|---|
| MAMP (235.30) | 0.136 | 0.897 | 0.123 | 0.051 |
| free MFT (114.16) | — | 1.85 | — | 0.096 |
| thiazole **6** (143.20) | 8.01 | 1.14 | 10.29 | 90.0 |
| methoxymethylpyrimidine (153.18) | 2.00 | — | not in Table 3 | not in Table 3 |
| MFT disulfide (226.31) | — | 0.376 | not in Table 3 | not in Table 3 |

Two checks that bound how far these numbers can be pushed:

- **Model II does not close a mass balance.** 211 mg/mol of free MFT is 1.85 mmol/mol × 3.32 mmol
  = 6.1 µmol recovered out of 8 760 µmol charged (0.07 %); MAMP took 0.897 × 3.32 = 3.0 µmol, i.e.
  **0.034 % of the MFT charged**. Either nearly all the MFT was lost to something unmeasured (volatility,
  oxidation, extraction) or the Table 2 figures are not absolute; no internal standard is stated for I/II.
  The tables are semi-quantitative GC-MS with no response factors, no replicates and no error bars.
- **In the thiamin-only aqueous pot (III), MFT bound in MAMP (0.123 mmol/mol) is reported where free MFT
  is not reported at all.** With cysteine (IV), free MFT 0.096 vs bound 0.051 mmol/mol. Thiazole release
  rose 8.7x with cysteine (12894 / 1474), so cysteine also changed the cleavage rate itself, by an
  unstated mechanism (pH was not reported).

### Mechanism (Fig. 2, p. 4057; text pp. 4057-4058) — proposed, not measured

Following Zoltewicz (bisulfite and amine substitution of thiamin): protonation of the pyrimidine,
OH⁻ attack, departure of the thiazole **6**, giving a stabilised pyrimidinyl-methylene cation
intermediate (**5**), which is captured by a nucleophile; with MFT (**7**) it gives **8** and then MAMP by
loss of water. Evidence for capture rather than rearrangement: adding MFT raised MAMP from 32 to 211
mg/mol (6.6x, Table 2), and in methanol the solvent-capture product (the methoxymethylpyrimidine, 307
mg/mol) appears without MFT and disappears with it. Cysteine, "a nucleophilic competitor of 7", lowered
MAMP (29 → 12) and raised MFT (blank → 11) (Table 3). The authors conclude that MFT formation from heated
thiamin is reduced because MFT is used up making MAMP (p. 4057). Odour: MAMP "had only a very weak sulfury odor" (p. 4057); no threshold given.

## 3. What it means for the model

1. **The paper establishes that an MFT sink by thiamin exists, qualitatively.** The electrophile is not
   intact thiamin but the pyrimidinyl-methyl fragment released when the thiazole leaves (the same bridge
   cleavage Mulley 1975 cites). MFT is trapped as a stable thioether that is nearly odourless. Confidence
   that the channel is real: high (isolated, structure-proven, synthesised, and responds to added MFT and
   to a competing nucleophile as the mechanism predicts). Confidence in any magnitude: low — the paper
   gives **no rate constant, no time course and no temperature dependence**, and its yields are
   semi-quantitative with an unclosed balance.
2. **The model has no such sink.** Read-only grep of `src/kinetic_core/sulfur.py` and
   `species_sulfur.py` (2026-10-09): no species or step for MAMP, the pyrimidine fragment, a
   furylthiomethyl thioether, or the thiazole **6**. Thiamine (`THI`, 12 C, 4 N) has exactly two fates,
   `r_thi_hmp` (→ HMP + FRAG_C 7 + FRAG_N 4) and `r_thi_mesh` (→ MeSH + FRAG_C 11 + FRAG_N 4); the
   pyrimidine half is discarded into inert fragments. The only MFT thioether channel,
   `ch_thioether_mft` (MFT + MELE → BND), drains the lumped matrix-electrophile pool `MELE`, whose only
   producers are the deoxyosone branches `ch_mele_from_dpo/tdp/ddp`; thiamine feeds no MELE.
3. **Where it would matter, in scaling terms (inference, not from the paper).** If capture follows the
   Fig. 2 mechanism, the electrophile is made at the thiamin cleavage rate, so the MFT drain per unit of
   MFT scales with the thiamin concentration. Jhoo's pots are 0.22-0.30 M thiamin. The panel pot the
   repo already notes (Zhang, 44.5 mmol/L thiamine; `ph_state.py`) is only 5x below Jhoo III, so the sink
   could be non-negligible there and in any thiamine-fortified model system. Beef at 0.0004-0.004 mM
   (`lombardiboccia2005_extraction.md`) is about 5 orders below, where it should be negligible. This also
   cuts against reading a thiamine-rich pot's MFT yield as the HMP route's yield: the model would
   attribute any MAMP-trapped MFT to a smaller `k_thi_hmp` or a larger MFT decay.
4. **The thiazole fate is missing too.** Intact thiazole **6** is the largest product in every model
   system (1.1-90 mmol/mol thiamin, Tables 2-3 derived); the model's thiamine has no thiazole-releasing branch that skips HMP.
   Mulley's total-loss rate is therefore a ceiling on `k_thi_hmp + k_thi_mesh`, and the HMP share of
   thiamine cleavage is unconstrained by either paper.

## What it does not give

- No rate constant for MAMP formation or MFT capture, no time course, no Ea; one temperature per system
  (80 °C for I-IV, oil bath 110 °C for the preparative run).
- No pH for I-IV; no pH measured after heating for any run.
- No detection limits, response factors, replicates or error bars; Table 2 has no stated internal
  standard; the Model II MFT balance does not close (0.07 % recovered).
- No quantity of free MFT in the preparative pH 6.5 buffer run, so no MAMP : MFT ratio under the one set
  of conditions that resembles a food pH.
- No HMP, no 3-mercapto-2-pentanone and no H₂S reported in any table.
- No odour threshold for MAMP ("very weak sulfury odor" only).
- No food matrix and nothing at food-level thiamin concentrations.
- The solvent of III/IV is stated two ways (water in the methods, methanol in the discussion).
