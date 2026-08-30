# Martins, Marcelis & van Boekel 2003 (10.1016/S0008-6215(03)00173-3) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/martins2003b.pdf` (13 pp.). Read method: both — `pdftotext -layout`
text layer for the body and Table 1, plus page rasters (PDF pp. 6, 7, 8, 9 = journal pp. 1656–1659)
read with the Read tool for all eight figures, which carry the bulk of the quantitative content.

---

## 0. PAPER IDENTITY — **MATCHES**, but the sibling keys around it DO NOT

> **IDENTITY, STATED IN BOLD AT THE TOP AS REQUIRED.**
> **`data/articles/martins2003b.pdf` IS PART I (Reaction mechanism), Carbohydr. Res. 338 (2003)
> 1651–1663.** This is correct as expected.
>
> **BUT THE WAVE BRIEF'S SIBLING ASSIGNMENTS WERE WRONG AND MUST BE RESTATED HERE:**
> - **`martins2003.pdf` is Part II (Kinetic analysis), Carbohydr. Res. 338 (2003) 1665–1678 — this
>   is where the rate-constant table lives.** It was *expected* to be `martins2003c.pdf`. It is not.
> - **`martins2003c.pdf` is a DIFFERENT paper entirely**: "Melanoidins extinction coefficient in
>   the glucose/glycine Maillard reaction", Food Chemistry 83 (2003) 135–142. It is neither Part I
>   nor Part II.
>
> A second, independent hazard: **the authors' own lettered self-citations are permuted relative
> to our repo keys, and permuted differently in different papers.** In Martins & Van Boekel (2005)
> Food Chem. **92**: "2003a" = Part II, "2003b" = extinction coefficient, "2003c" = Part I. In
> Martins & Van Boekel (2005) Food Chem. **90**: "2003a" = Part I, "2003b" = Part II,
> "2003c" = extinction coefficient. **Never resolve a Martins 2003 citation by its letter suffix.
> Resolve it by journal + volume + page.**

| field | value |
|---|---|
| authors | Sara I.F.S. Martins,ᵃ Antonius T.M. Marcelis,ᵇ Martinus A.J.S. van Boekelᵃ,* |
| title | "Kinetic modelling of Amadori N-(1-deoxy-D-fructos-1-yl)-glycine degradation pathways. Part I — Reaction mechanism" |
| venue | Carbohydrate Research |
| volume/pages/year | 338 (2003) 1651–1663 |
| DOI | 10.1016/S0008-6215(03)00173-3 |
| received / accepted | Received 11 December 2002; accepted 11 April 2003 |
| affiliation | ᵃ Product Design and Quality Management Group, Dept. of Agrotechnology and Food Sciences, Wageningen University, P.O. Box 8129, 6700 EV Wageningen, The Netherlands; ᵇ Organic Chemistry Group, same department and address |
| PDF character | born-digital (Elsevier / Acrobat Distiller 4.05 for Windows); clean text layer, but "±" renders as `9/` and "−"/"–" as `/`; **all eight figures are image objects with no extractable data** |
| corresponding author | tiny.vanboekel@wur.nl |
| funding | Fundação para a Ciência e Tecnologia (Portugal) |
| notable acknowledgement | DFG material "kindly offered by Drs. Blank and Davidek, from Nestlé Research Center Lausanne" |

### PRE-EXISTING-COVERAGE FLAG — DOUBLE-COUNTING HAZARD

`docs/reference/FIT_HOLDOUT_DECLARATION.md` (read-only; Amendment 1, 2026-08-28, lines ~209–215)
declares roles for "Martins 2003 Wageningen thesis (edepot.wur.nl/121418)" Tables 4.2.3, 4.1.1,
3.3.1 and Ch. 6. **These four journal papers are the thesis chapters.**

**My mapping proposal for this paper (proposal, NOT a verified fact — I have not opened the
thesis PDF):**

| declared thesis item | my proposed mapping | confidence | verdict |
|---|---|---|---|
| Table **4.1.1** "glycine-release yields (65–95%)", role **FIT** | **This paper's Table 1** (p. 1657), "Glycine yield from DFG heated in phosphate buffer (0.1 M) at different reaction conditions". The declared range **65–95%** is a near-exact description of the pH 6.8 columns B and D of that table (B spans 71.02–93.85; D spans 83.92–95.82). | **HIGH** | **DIRECT DOUBLE-COUNT. Same measurements, re-published.** If Table 4.1.1 is already FIT, then this paper's Table 1 is already FIT and must not be re-declared. |
| Table **4.2.3** "reverse rates (Amadori → parent sugar)" | not in this paper — this paper has **no rate constants at all** | high | n/a here; but see the blocking note in `martins2003_extraction.md` §4.4 — I could not find such reverse rates in *any* of the four papers, nor in `martins2005.pdf` |
| Table 3.3.1 melanoidin ε | maps to `martins2003c.pdf` Table 1 | high | n/a here (this paper only *cites* ε = 0.64 ± 0.03) |
| Ch. 6 pH ladder | maps to `martins2005b.pdf` Tables 2 and 3 | high | n/a here |

**Second, intra-wave double-count:** **Part II (`martins2003.pdf`) Figs. 2 and 4–7 replot the same
experimental points as this paper's Figs. 1–8**, with model curves overlaid. Part I and Part II
are **one experiment reported twice.** Part II's Table 3 is a regression *on* these data, so it is
not independent of them either. Assign the time-course role **once**, here, and cross-reference.
Assigning it in both dossiers would be exactly the rule-1 disjointness violation that Amendment 2
of the declaration had to retro-fix for Cerny 2007.

---

## 1. ONE-PARAGRAPH VERDICT

This is the experimental half of the Martins DFG-degradation pair. It gives the model a complete
2 × 2 (100/120 °C × initial pH 5.5/6.8) time-resolved dataset for the thermal decomposition of the
**isolated Amadori compound** with ten simultaneously quantified responses (DFG, glycine, 3-DG,
1-DG, methylglyoxal, glucose, mannose, fructose, formic acid, acetic acid, melanoidins) plus a
measured in-run pH drop. **It contains NO rate constants, NO activation energies, and NO fitted
parameters of any kind** — the word "kinetic" in the title refers to the *forthcoming* analysis in
Part II. Only **one table** is printed (Table 1, glycine yields); **every other number in the paper
must be digitised from a figure.** There is **no mass-balance / carbon-balance figure in this
paper** — the mass balance the wave brief asked for lives in **Part II Fig. 3** (100 and 120 °C at
pH 6.8 only) and is a *molar* balance against [DFG]₀, not a carbon balance; it is transcribed in
`martins2003_extraction.md` §8.3. What this paper *does* uniquely give is (a) the glycine-release
yield table, (b) the proposed reaction scheme as an explicit primary/secondary route list, and
(c) two auxiliary control experiments (Glc alone and Fru alone, 120 °C pH 6.8) whose results are
reported **qualitatively only, "results not shown"** — an important negative.

---

## 2. SYSTEM DEFINITION — verbatim

| variable | value as printed | anchor |
|---|---|---|
| reactant | N-(1-deoxy-D-fructos-1-yl)-glycine (DFG), compound **1** | §2.2, p. 1653 |
| DFG synthesis | Glycine (5.25 g, 0.07 mol), D-Glc (50.3 g, 0.28 mol), Na₂S₂O₅ (7.6 g, 0.04 mol) in 1:1 MeOH–water (80 mL), stirred 20 min at rt, refluxed 8 h; Dowex 50W×8 (H⁺) ion exchange, eluted with 0.1 mol/L ammonia; "The final yield was **19%** based on the amount of glycine used." | §2.2, p. 1653 |
| DFG purity | "Anal. Calcd for C₈H₁₅NO₇: C, 40.50; H, 6.37; N, 5.91. Found: C, **40.13**; H, **6.4**; N, **5.89**. mp: **160–163 °C**." Single peak on the amino acid analyser and single peak by HPLC. Additional DFG donated by Nestlé Research Center Lausanne. | §2.2, p. 1653 |
| initial concentration | "Compound 1 (0.237 g, 10 mmol) was dissolved in phosphate buffer (100 mL …)" → **10 mmol/l** | §2.4, p. 1654 |
| buffer | "phosphate buffer (100 mL, **0.1 M** K₂HPO₄/KH₂PO₄, pH 6.8 and 5.5)" | §2.4, p. 1654 |
| pH | initial 5.5 and 6.8. **NOT hot-pH corrected** and **NOT held constant** — the paper reports the drift explicitly (§3.1.6 and Fig. 7): at initial pH 6.8, drop = **0.21** after 120 min at 120 °C and **0.19** after 180 min at 100 °C; at initial pH 5.5, "Due to the low buffer capacity … the pH drop for the same reaction conditions was **0.8 and 0.9**, respectively" | §3.1.6, pp. 1658–1659 |
| temperature | 100 °C "for a maximum of 4 h" and 120 °C "for a maximum of 3 h", **oil bath** | §2.4, p. 1654 |
| vessel | "screw-capped glass tubes (Schott, 16 × 160 mm)" | §2.4, p. 1654 |
| headspace / atmosphere | **not stated.** No inerting, no headspace volume, no purge reported. | — |
| volume | 100 mL prepared; sub-samples withdrawn at predetermined times | §2.4, p. 1654 |
| agitation | **not stated** | — |
| quench | "At predetermined heating times, samples were taken and immediately cooled in **ice water**, prior to analyses." | §2.4, p. 1654 |
| replication (n) | "Each reaction mixture was prepared, heated and analysed in **at least duplicate**." | §2.4, p. 1654 |
| error-bar definition | **NEVER DEFINED IN THIS PAPER.** Error bars appear on Figs. 1, 2, 3, 4, 5, 6, 8 but neither the captions nor §2 nor §3 say whether they are SD, SEM, or range. **This is a real gap: unstated.** (Part II §3.4 only says "the scatter in replicates was not very high".) | — |

### Analytical methods, verbatim details

**DFG (1) and glycine — §2.5, p. 1654.** Reaction mixtures diluted 1:10 with water, pH adjusted to
2.2 with loading buffer (LB, Biochrom) in the proportion 100 µL diluted sample : 900 µL LB.
Loaded: 20 µL sample + 20 µL internal standard (**glutamic acid; printed as "25 mmol/mL"**, which
is physically impossible and is a **printing/unit error — most plausibly 25 nmol/µL, i.e. 25 mM**)
+ 20 µL LB. Biochrom 20 Amino Acid Analyser (Pharmacia, England). Separation program, verbatim:
lithium citrate pH 2.8 / 16 min / 80 °C → Li citrate pH 3.0 / 21 min / 75 °C → Li citrate pH 3.15 /
12 min / 75 °C → Li citrate pH 3.55 / 21 min / 75 °C → lithium hydroxide / 6 min / 85 °C → Li
citrate pH 2.8 / 8 min / 85 °C → Li citrate pH 2.8 / 20 min / 75 °C → Li citrate pH 2.8 / 5 min /
80 °C. Post-column ninhydrin derivatisation, photometric detection at **570 nm**.
t_R: DFG **9.5 min**; glycine **38.8 min**. "Quantification was done using **external standards**."

**1-DG (5), 3-DG (3), methylglyoxal (8) — §2.6, p. 1654.** 1 mL sample + 1 mL water + 2 mL
methanolic 1,2-diaminobenzene (**1 mol/L**), held **overnight at rt (25 °C)**; "With previous
experimental tests it was established that this treatment ensured a **complete derivatization** of
the target compounds into their quinoxaline form." RP-HPLC, gradient from 10:90 (v/v)
acetonitrile : ammonium acetate buffer (**pH 3.5, 20 mmol/L**), MeCN raised to **30% within 50 min**,
flow **0.4 mL/min**, detection **320 nm**, column temperature raised to **40 °C** to optimise
separation. External standards: 2-methyl-3-(1,2,3-trihydroxypropyl)quinoxaline (5′) and
2-methylquinoxaline (10′); "**Both quinoxalines have the same extinction coefficient.**"
t_R: 5′ **11.3 min**; 2-(2,3,4-trihydroxybutyl)quinoxaline (3′) **12.4 min**;
2-methylquinoxaline **38.8 min**.
The 5′ reference standard was synthesised in-house (§2.3) in **0.02% yield**, ¹H NMR (200 MHz,
CD₃OD) δ 2.86 (s, 3H, –CH₃), 3.81–3.92 (m, 2H, –CH₂OH), 4.08 (m, 1H, –CH(OH)–), 5.10 (d, 1H,
–CCH(OH)–), 7.77 (m, 2H, –CH–), 7.98 (m, 1H, –CH–), 8.096 (m, 1H, –CH–).

**Glc (10), Man (11), formic acid, acetic acid — §2.7, p. 1654.** 1:10 dilution; HPLC on ION-300
ion-exchange (Interaction Chromatography, San Jose CA); eluent **2.5 mmol/L H₂SO₄** in Millipore
water, flow **0.4 mL/min**, column **85 °C**. Sugars by refractive index: Glc **15.1 min**,
Man **16.4 min**. Acids by A₂₁₀: formic **21.6 min**, acetic **23.7 min**. All by external standards.
**CRITICAL LIMITATION, stated by the authors (§3.1.5): on this column "Man and Fru have the same
retention time"** — so every Man value from the ION-300 method is really Man + Fru.

**Glc / Man / Fru resolved — §2.8, p. 1655.** 1:100 dilution; Dionex (Sunnyvale CA) CarboPac PA100;
gradient 16/0/84 (v/v) 0.1 M NaOH / 1 M NaOAc in 0.1 M NaOH / water for 20 min → NaOAc to 0/100/0
isocratic 10 min → 100/0/0 in 1 min, isocratic 5 min → back to start, 25 min. Flow **1 mL/min**.
ED-40 electrochemical detector. t_R: Glc **12.95 min**, Man **13.9 min**, Fru **16.59 min**.
External standards.

**Heterocycles (HMF, HHMF, DDMP) — §2.9, p. 1655.** RP-HPLC (Lichrosorb RP-18, Merck), eluent
**7.5:100 MeOH–water**, flow **0.6 mL/min**, UV **280 nm**. "Only HMF was quantified by external
calibration line."

**Melanoidins — §2.10, p. 1655.** A₄₇₀ (Pharmacia Biotech, Uppsala), diluted with demi water when
necessary, converted by Lambert–Beer using
**ε = 0.64 ± 0.03 l/mmol·cm at 470 nm** [C, cited to ref. 17 = `martins2003c.pdf`].
"The concentration of melanoidins is thus expressed as **moles of sugar incorporated in the brown
polymers**."

---

## 3. TABLE 1 — GLYCINE YIELD (the only table in the paper) — complete transcription

**Anchor: Table 1, p. 1657 (PDF p. 7).** Title as printed:
> "Glycine yield from DFG heated in phosphate buffer (0.1 M) at different reaction conditions"

Column headers as printed: `Heating time (min)` | `Glycine yield (%) ᵃ` spanning four sub-columns
`A (T 100 °C pH 5.5)` | `B (T 100 °C pH 6.8)` | `C (T 120 °C pH 5.5)` | `D (T 120 °C pH 6.8)`

Footnotes as printed:
> "n.a.: not analyzed."
> "ᵃ The yield of glycine was calculated based on the amounts of DFG reacted."

**All values [M] — measured (a ratio of two measured concentrations; the authors report it as data,
not as a fit). Units: % (mol glycine formed per mol DFG reacted × 100).**

| Heating time (min) | A (T 100 °C pH 5.5) | B (T 100 °C pH 6.8) | C (T 120 °C pH 5.5) | D (T 120 °C pH 6.8) |
|---|---|---|---|---|
| 5 | 43.40 | n.a. | 25.87 | n.a. |
| 10 | n.a. | n.a. | 62.05 | 83.92 |
| 15 | 57.89 | 71.02 | 59.56 | 94.72 |
| 30 | 53.61 | 93.85 | 69.79 | 95.82 |
| 45 | 78.23 | 87.44 | 66.20 | 92.36 |
| 60 | 79.69 | 84.91 | 71.95 | 84.41 |
| 75 | 81.46 | n.a. | n.a. | n.a. |
| 90 | 72.39 | 78.42 | 74.20 | 84.05 |
| 120 | 65.03 | 76.59 | 70.18 | 86.39 |
| 150 | n.a. | 80.14 | n.a. | n.a. |
| 180 | 78.31 | 81.45 | n.a. | n.a. |

**No error column is printed.** No SD, no n per cell. The 2-decimal precision is not supported by
any stated uncertainty and should not be read as 4-significant-figure accuracy.

### Derived summary [Z] — my computation from the cells above, never printed by the authors

| quantity [Z] | A (100/5.5) | B (100/6.8) | C (120/5.5) | D (120/6.8) |
|---|---|---|---|---|
| min yield (%) | 43.40 (5 min) | 71.02 (15 min) | 25.87 (5 min) | 83.92 (10 min) |
| max yield (%) | 81.46 (75 min) | 93.85 (30 min) | 74.20 (90 min) | 95.82 (30 min) |
| yield at final tabulated time (%) | 78.31 (180 min) | 81.45 (180 min) | 70.18 (120 min) | 86.39 (120 min) |
| mean over tabulated times (%) | 67.8 | 82.0 | 62.6 | 89.5 |

The declaration's "**65–95%**" descriptor for thesis Table 4.1.1 is consistent with the pH 6.8
columns (B: 71.0–93.9; D: 83.9–95.8) and with the pH 5.5 columns' *upper* portion, but **not**
with the two low early points (25.87 and 43.40). If the declaration's range was meant to be
inclusive, it under-reports the spread.

### Authors' interpretation of Table 1 — verbatim, [M]

> "at the early stage of the reaction the yield of glycine increased with time … However, as the
> reaction proceeded a decrease was observed. Moreover, in all the studied systems the decrease in
> the yield of glycine was followed by an increase again. Also, the yield of glycine increased with
> increasing pH, independently of the temperature." (§3.1.3, p. 1656)

Mechanistic reading (§3.2.2, p. 1660):
> "In theory, **1 mol of DFG should yield 1 mol of glycine**. … This indicates that glycine is
> first liberated and then reacts with other compounds present in the reaction mixture (e.g.,
> α-dicarbonyls and hydroxycarbonyls). These reactions may include: (i) the formed α-dicarbonyls
> can further react with glycine through **Strecker degradation** …; (ii) formation of
> **pyrazinones and pyrazines**; (iii) **chain elongation** of reactive C-2 and C-3 α-dicarbonyls
> by one carbon unit originating from C-2 atom of glycine; (iv) **incorporation of glycine into
> melanoidins**."

Abstract-level headline yields, verbatim from the Abstract, p. 1651 [M]:
> "Together with glycine, acetic acid was the main end product formed. Values of **83 and 55 mol%**
> were obtained, respectively."
**The conditions for these two headline numbers are NOT stated in the abstract, and I cannot
locate a body-text or table cell that reproduces "83" and "55" exactly.** 83 mol% is close to the
pH 6.8 glycine plateau (B ≈ 81.5, D ≈ 86.4). 55 mol% has no matching cell anywhere; from Fig. 6
the acetic acid at 100 °C/pH 6.8 reaches ≈4.8 mmol/l on 10 mmol/l DFG ≈ 48 mol%, and at
120 °C/pH 6.8 ≈5.6 mmol/l ≈ 56 mol% [Z from fig]. **Best reconstruction: 83 mol% glycine and
55 mol% acetic acid both at pH 6.8, temperature ambiguous. Flagged as an under-specified abstract
claim; do not quote "83 and 55 mol%" without stating the conditions are not given.**

---

## 4. THE PROPOSED REACTION SCHEME, STEP BY STEP

### 4.1 Scheme 1 (p. 1652) — the mechanistic chemistry (numbered compounds)

**Anchor: Scheme 1, p. 1652 (PDF p. 2).** Caption verbatim: "Amadori compound
N-(1-deoxy-D-fructos-1-yl)glycine degradation pathways: enolization and retro-aldolization."
The scheme is a structure drawing; the *pathway* content is given in the body text (§1, pp.
1651–1653) and is transcribed here as arrows. Numbered compounds used throughout the paper:

| no. | compound |
|---|---|
| **1** | N-(1-deoxy-D-fructos-1-yl)-glycine (DFG), the Amadori compound |
| **2** | 1,2-enaminol |
| **3** | 3-deoxy-2-hexosulose (3-DG) |
| **4** | 2,3-enaminol |
| **5** | 1-deoxy-2,3-hexodiulose (1-DG) |
| **6** | 1-glycine-1-deoxy-D-glyceraldehyde |
| **7** | glyceraldehyde |
| **8** | methylglyoxal (MG) |
| **9 / 9′** | Schiff base (two forms) |
| **10** | D-glucose |
| **11** | D-mannose |
| **12** | 2,3-endiol |
| **13** | D-fructose |
| **3′** | 2-(2,3,4-trihydroxybutyl)quinoxaline (derivative of 3) |
| **5′** | 2-methyl-3-(1,2,3-trihydroxypropyl)quinoxaline (derivative of 5) |
| **8′ / 10′** | 2-methylquinoxaline (derivative of 8) |

Pathways as lettered in Scheme 1, verbatim from §1:
- **Pathway A — 1,2-enolization:** 1 → 2 (1,2-enaminol) → **3** (3-DG), "accompanied by amino acid
  release".
- **Pathway B — 2,3-enolization:** 1 → 4 (2,3-enaminol) → **5** (1-DG), likewise with amino acid
  release.
- **Pathway C — retro-aldol:** "1 generates 1-glycine-1-deoxy-D-glyceraldehyde (**6**) and
  glyceraldehyde (**7**) through a retro-aldol cleavage at **C-3 to C-4**. Compound 7 can react with
  free glycine and produce more of compound 6 that subsequently undergoes a β-elimination to form
  methylglyoxal (**8**) and release glycine. Moreover, 8 can also be produced from 7 through the
  catalytic action of the amino acid."
- **Pathway D — reversal (the reversibility hypothesis):** "the Amadori rearrangement is expected
  to be a reversible process, that is to undergo a non-enzymatic reversal into enolamines and
  subsequently into free sugars and amino acid, as shown in pathway D."
  **This pathway is DRAWN but never quantified. No rate constant is attached to it in this paper.**
- Side branches noted: transition-metal-catalysed oxidation of the 1,2-enaminol (2) to glucosones;
  elimination of the C-4 OH of the 2,3-enaminol (4) to give 1-amino-1,4-dideoxy-2,3-diulose.

### 4.2 Scheme 3 (p. 1661) — organic-acid formation
**Anchor: Scheme 3, p. 1661 (PDF p. 11).** Caption verbatim: "Reaction scheme for organic acids
formation from primary thermal degradation products of DFG (**adapted from Ginz and co-workers³⁴**)."
**[C] — this scheme is adapted from Ginz, Balzer, Bradbury & Maier, Eur. Food Res. Technol. 2000,
211, 404–410, not derived here.** Its content per the body text: Glc, Man and Fru together with
1-DG and 3-DG are aliphatic-acid precursors; formic acid from **C-1 to C-2 cleavage of 3-DG**,
acetic acid from **C-2 to C-3 cleavage of 1-DG** or from cleavage of triose intermediates such as
MG (§1, p. 1653).

### 4.3 §3.3 — the authors' final PRIMARY / SECONDARY route list, verbatim

This is the deliverable the wave brief asked for as "the proposed reaction scheme step by step."
**Anchor: §3.3.1 and §3.3.2, pp. 1661–1662.** Reproduced verbatim, bullet for bullet:

**3.3.1. Primary routes.**
1. "DFG degradation through enolization occurs through two intermediates, designed as E₁ and E₂,
   which can be the Schiff's base, the cation form of the Schiff's base, the 1,2-enaminol or the
   2,3-enaminol. These intermediates have not been isolated yet from the Maillard reaction due to
   their reactivity, however, according to previous studies ARPs undergo 1,2- or 2,3-enolizations
   while the amino acid moiety is still attached. **E₁ is favoured at lower pH whereas E₂ is
   favoured at higher pH.**"
2. "The intermediates E₁ and E₂ by release of the amino acid lead to the formation of **3-DG and
   1-DG**, respectively."
3. "**Glc, Man and 3-DG are formed through the same intermediate (E₁) whereas Fru is formed by
   E₂.** Moreover, Man and Glc can isomerise into each other as well as degrade into Cₙ (n ≤ 6)
   carbonyl compounds."
4. "DFG degradation through **retro-aldol cleavage leads to MG formation, with amino acid release**."
5. "DFG, 3-DG and 1-DG due to their reactive functional group can easily degrade to produce
   reactive Cₙ carbonyl compounds."
6. "**Formic acid** and **acetic acid** are formed from DFG degradation pathways **through E₁ by
   3-DG dicarbonyl cleavage** and **through E₂ by 1-DG dicarbonyl cleavage**, respectively."
7. "**Melanoidins formation results from the interaction of Cₙ carbonyl compounds with glycine.**"

**3.3.2. Secondary routes.**
1. "Products formed by cyclization/condensation of the deoxyosones that include **HMF, HHMF and
   DDMP**."
2. "Incorporation of glycine in **pyrazinones and pyrazines**, as well as in chain elongation of
   α-dicarbonyls. Also degradation of glycine through **Strecker degradation** to produce the
   corresponding amines, carboxylic acids and Strecker aldehydes, which can be included in the
   melanoidins formation."
3. "**MG formation from deoxyosones retro-aldol reaction** and its degradation into carboxylic acids."
4. "**Isomerisation of Glc and Man into Fru**, as well as **direct degradation of Glc into 3-DG**."
5. "**Fru formation from E₁** as well as its degradation into reactive Cₙ carbonyl compounds."
6. "**Carboxylic acids formation from sugars enolization** or scission products from Cₙ carbonyl
   compounds."

**Note on continuity:** Part II's Scheme 1 is this list drawn as a graph; Part II's Scheme 2 (M1)
and Scheme 3 (M2) are its fittable simplifications. **Secondary route 4 ("direct degradation of Glc
into 3-DG") is precisely the step that Part II promoted from secondary to primary as step 13 when
M1 was replaced by M2.** That promotion is traceable and is a genuine model-structure finding.

---

## 5. FIGURE 1 — DFG degradation

**Anchor: Fig. 1, p. 1656 (PDF p. 6).** Caption verbatim: "DFG thermal degradation at 100 (△,
pH 5.5; ▲, pH 6.8) and at 120 °C (○, pH 5.5; ●, pH 6.8)."
Axes as printed: y = `DFG mmol/l`, 0–12, ticks every 2; x = `Time (min)`, 0–200, ticks every 50.
Error bars present, **definition unstated**.

Digitised from the PDF p. 6 raster. **All values [fig].** Read precision ≈ ±0.3 mmol/l
(≈ ±2.5% of full scale); marker overlap near zero limits resolution below ~0.5 mmol/l.

| time (min) | △ 100 °C pH 5.5 | ▲ 100 °C pH 6.8 | ○ 120 °C pH 5.5 | ● 120 °C pH 6.8 |
|---|---|---|---|---|
| 0 | ~9.6 | ~9.6 | ~9.8 | ~9.8 |
| 5 | ~9.5 | ~8.5 | ~8.7 | ~4.6 |
| 10 | ~9.3 | ~7.4 | ~7.5 | ~2.0 |
| 15 | ~9.2 | ~6.4 | ~6.6 | ~1.0 |
| 30 | ~8.7 | ~4.4 | ~4.7 | ~0.1 |
| 45 | ~8.2 | ~3.0 | ~3.2 | ~0 |
| 60 | ~7.7 | ~1.9 | ~2.2 | ~0 |
| 90 | ~6.7 | ~0.7 | ~1.0 | ~0 |
| 120 | ~5.9 | ~0.2 | ~0.4 | ~0 |
| 150 | ~5.1 | unreadable (overlapping at ~0) | unreadable | ~0 |
| 180 | ~4.3 | ~0 | ~0 | — |

Authors' own quantitative statements [M], which should be preferred over my [fig] reads where they
overlap (§3.1.1, p. 1655):
> "After heating for **30 min at 120 °C and pH 6.8, 1 was completely degraded** whereas after
> **180 min at 100 °C and pH 5.5, only approximately 57% of the initial concentration had reacted**.
> Moreover, when comparing the system at 100 °C and pH 6.8 with the system at 120 °C and pH 5.5, a
> **similar degradation rate** was observed. These results indicate that **decreasing the initial
> reaction pH with 1.3 units or increasing the temperature with 20 °C has the same effect on the
> DFG degradation.** However, this observation is only valid in the pH range between 5 and 7."

My [fig] read gives 9.6 → 4.3 at 100 °C/pH 5.5 over 180 min = **55% reacted**, against the authors'
"approximately 57%" — agreement within my stated read precision. Good cross-check on the
digitisation.

---

## 6. FIGURE 2 — glycine formation

**Anchor: Fig. 2, p. 1656 (PDF p. 6).** Caption verbatim: "Glycine formation during thermal
degradation of DFG at 100 (△, pH 5.5; ▲, pH 6.8) and 120 °C (○, pH 5.5; ●, pH 6.8)."
Axes as printed: y = `Glycine mmol/l`, 0–10, ticks every 2; x = `Time (min)`, 0–200, ticks every 50.
Error bars present, definition unstated.

Digitised from PDF p. 6 raster. **All [fig].** Read precision ≈ ±0.3 mmol/l.

| time (min) | △ 100 °C pH 5.5 | ▲ 100 °C pH 6.8 | ○ 120 °C pH 5.5 | ● 120 °C pH 6.8 |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 5 | ~0.3 | ~1.4 | ~1.2 | ~4.5 |
| 10 | ~0.5 | ~2.5 | ~2.0 | ~6.0 |
| 15 | ~0.7 | ~3.4 | ~2.6 | ~6.8 |
| 30 | ~1.2 | ~5.4 | ~4.0 | ~7.6 |
| 45 | ~1.9 | ~6.4 | ~5.0 | ~7.7 |
| 60 | ~2.4 | ~7.0 | ~5.9 | ~7.8 |
| 90 | ~3.2 | ~7.5 | ~6.8 | ~8.0 |
| 120 | ~3.8 | ~7.7 | ~7.1 | ~8.2 |
| 150 | ~4.3 | ~7.9 | unreadable | — |
| 180 | ~4.7 | ~8.1 | ~7.3 | — |

Consistency check [Z]: at 120 °C / pH 6.8, [Gly] ≈ 8.2 mmol/l against [DFG]₀ ≈ 9.8 mmol/l = 84 mol%,
matching Table 1 column D (84.05–86.39% at 90–120 min) and the abstract's "83 mol%". **The
digitisation is self-consistent with the table.**

---

## 7. FIGURE 3 — deoxyosone formation

**Anchor: Fig. 3, p. 1657 (PDF p. 7).** Caption verbatim: "Deoxyosones formation during thermal
degradation of DFG at 100 (A) and 120 °C (B): 1-DG (△, pH 5.5; ▲, pH 6.8) and 3-DG (○, pH 5.5; ●,
pH 6.8)."
Axes as printed: **panel A** y = `Deoxyosones mmol/l`, 0–0.5, ticks 0.05; x = `Time (min)` 0–200.
**Panel B** y = `Deoxyosones mmol/l`, 0–0.5, ticks 0.05; x = `Time (min)` 0–150.
Error bars present, definition unstated.

Digitised from PDF p. 7 raster. **All [fig].** Read precision ≈ ±0.015 mmol/l — this is a
**low-signal figure** and the four series overlap heavily below 0.1 mmol/l.

**Panel A — 100 °C**

| time (min) | △ 1-DG pH 5.5 | ▲ 1-DG pH 6.8 | ○ 3-DG pH 5.5 | ● 3-DG pH 6.8 |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 5 | ~0.05 | ~0.05 | ~0.06 | ~0.06 |
| 15 | ~0.07 | ~0.045 | ~0.13 | ~0.04 |
| 30 | ~0.09 | ~0.035 | ~0.20 | ~0.025 |
| 45 | ~0.085 | ~0.025 | ~0.26 | ~0.02 |
| 60 | ~0.08 | ~0.02 | ~0.31 | ~0.015 |
| 90 | ~0.07 | ~0.01 | ~0.36 | ~0.01 |
| 120 | ~0.06 | ~0.01 | ~0.40 | ~0.005 |
| 150 | ~0.05 | ~0.005 | ~0.43 | ~0.005 |
| 180 | ~0.04 | ~0.005 | ~0.47 | ~0.005 |

**Panel B — 120 °C**

| time (min) | △ 1-DG pH 5.5 | ▲ 1-DG pH 6.8 | ○ 3-DG pH 5.5 | ● 3-DG pH 6.8 |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 5 | ~0.06 | ~0.05 | ~0.13 | ~0.045 |
| 10 | ~0.07 | ~0.035 | ~0.28 | ~0.02 |
| 15 | ~0.08 | ~0.02 | ~0.37 | ~0.01 |
| 30 | ~0.055 | ~0.01 | ~0.43 | ~0.005 |
| 45 | ~0.03 | ~0.005 | ~0.40 | ~0 |
| 60 | ~0.02 | ~0.005 | ~0.33 | ~0 |
| 90 | ~0.01 | ~0 | ~0.24 | ~0 |
| 120 | ~0.005 | ~0 | ~0.16 | ~0 |

Authors' statements [M] (§3.1.4, p. 1656):
> "a quite high concentration of 1-DG and 3-DG was already present **after 5 min at pH 5.5 and
> temperature 100 °C (34%)**." *(The "34%" is not defined against a stated denominator — it reads as
> % of DFG reacted at that time. Ambiguous as printed.)*
> "Independently of the temperature, at lower pH 3-DG was present in larger amounts relatively to
> 1-DG. At pH 6.8 that difference decreased. An increase in pH seems to favour the 1-DG formation."
> "as the reaction proceeded the amount of **1-DG decreased much faster than the one of 3-DG**."
> "(DFG pKₐ is **8.2**)" — [C], quoted as literature context for 2,3-enolisation above the pKₐ.

And (§3.2.1, p. 1659):
> "**In the present study at 100 °C and pH 6.8 a ratio of 3.2 was found**" (3-DG : 1-DG), compared
> against "Hofmann and co-workers … at pH 7 found 3-DG in **four times** higher amounts than 1-DG"
> [C].

**INTERNAL INCONSISTENCY, FLAGGED:** my [fig] read of panel A at 100 °C / pH 6.8 gives
3-DG (●) ≈ 0.06 and 1-DG (▲) ≈ 0.05 at 5 min → ratio ≈ **1.2**, and the ratio *falls* below 1 after
~15 min. I cannot reproduce the stated **3.2** from the figure as captioned at any time point.
Either (a) the "3.2" refers to an integrated or early-time quantity not plotted, or (b) the ○/●
marker assignment in the caption is transposed. **The signal here is at the very bottom of the
figure's resolution, so I flag this rather than resolve it. Any wave using 3-DG:1-DG ratios at
pH 6.8 should treat the value as uncertain between ~1.2 [fig] and 3.2 [M].**

**Also literature-cited for context, [C]** (§3.2.1, p. 1659): heating 1-deoxy-1-propylamino-D-
fructose in phosphate buffer for 10 h under reflux gave 1-DG/3-DG = **8:5 at pH 4.5** and
**20:1 at pH 7** (Beck, Ledl & Severin, Carbohydr. Chem. 1988, 177, 240–243, ref. 24). Note this
*contradicts* the present study's 3-DG > 1-DG at every condition; the authors say so:
"In literature the results are however contradictory."

---

## 8. FIGURE 4 — methylglyoxal formation

**Anchor: Fig. 4, p. 1657 (PDF p. 7).** Caption verbatim: "Methylglyoxal formation during thermal
degradation of DFG at 100 (△, pH 5.5; ▲, pH 6.8) and 120 °C (○, pH 5.5; ●, pH 6.8)."
Axes as printed: y = `Methylglyoxal mmol/l`, 0–3.5, ticks 0.5; x = `Time (min)`, 0–140, ticks 20.
Error bars present, definition unstated.

Digitised from PDF p. 7 raster. **All [fig].** Read precision ≈ ±0.08 mmol/l.

| time (min) | △ 100 °C pH 5.5 | ▲ 100 °C pH 6.8 | ○ 120 °C pH 5.5 | ● 120 °C pH 6.8 |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 5 | ~0.10 | ~0.35 | ~0.75 | ~1.00 |
| 10 | ~0.18 | ~0.55 | ~1.20 | ~1.45 |
| 20 | ~0.30 | ~0.80 | ~1.70 | ~1.95 |
| 30 | ~0.42 | ~1.00 | ~2.05 | ~2.35 |
| 45 | ~0.60 | ~1.30 | ~2.40 | ~2.70 |
| 60 | ~0.72 | ~1.50 | ~2.60 | ~2.90 |
| 120 | ~1.15 | ~2.00 | ~2.95 | ~3.20 |

Authors' statements [M] (§3.1.4b, p. 1656):
> "Contrary to the temperature increase, the **increase of pH had almost no influence on the MG
> formation**. … However, the influence of temperature is quite significant. **An increase from 100
> to 120 °C more than doubled the formation of MG. No decrease in its amount was observed with
> time.**"

Cross-check [Z]: at 120 min, 100 °C pH 5.5 → 120 °C pH 5.5 is 1.15 → 2.95 = **2.6×**; at pH 6.8,
2.00 → 3.20 = **1.6×**. "More than doubled" holds at pH 5.5 but **not** at pH 6.8 on my read.
Flagged as a mild overstatement in the authors' text.

---

## 9. FIGURE 5 — parent sugar formation (Glc / Man / Fru)

**Anchor: Fig. 5, p. 1658 (PDF p. 8).** Caption verbatim: "Sugars formation during thermal
degradation of DFG at pH 5.5 (---) and pH 6.8 (···). **A and B heated at 100 °C; C and D heated at
120 °C.** D-Glc (●); D-Man (▲); D-Fru (×)."

**CAPTION AMBIGUITY, FLAGGED:** the caption names two line styles for two pHs and four panels for
two temperatures, but does not explicitly say which panel is which pH. **My reading, from the body
text plus panel content:** A = 100 °C pH 5.5, B = 100 °C pH 6.8, C = 120 °C pH 5.5,
D = 120 °C pH 6.8. Evidence: panels A and C contain **no × (Fru) series at all**, and the body text
says Fru was "observed in addition to Glc and Man, in samples heated at pH initial value 6.8"; and
the body says "at pH 5.5, Man was formed in higher amounts than Glc (Fig. 5 A and C)", which
matches A and C having ▲ above ●. **This is a reasoned assignment, not a printed one.**

Axes as printed: A y = `Sugars mmol/l` 0–0.5 (ticks 0.1), x 0–200; B y 0–0.7 (ticks 0.1), x 0–200;
C y 0–0.7 (ticks 0.1), x 0–150; D y 0–1 (ticks 0.2), x 0–150.
Error bars present on all panels, definition unstated.

Digitised from PDF p. 8 raster. **All [fig].** Read precision ≈ ±0.02 mmol/l.

**Panel A — 100 °C, pH 5.5** (no Fru detected)

| time (min) | ● Glc | ▲ Man | × Fru |
|---|---|---|---|
| 0 | 0 | 0 | not detected |
| 15 | ~0.02 | ~0.03 | — |
| 30 | ~0.06 | ~0.10 | — |
| 45 | ~0.10 | ~0.17 | — |
| 60 | ~0.14 | ~0.24 | — |
| 75 | ~0.18 | ~0.29 | — |
| 90 | ~0.21 | ~0.33 | — |
| 120 | ~0.26 | ~0.39 | — |
| 150 | ~0.30 | ~0.42 | — |
| 180 | ~0.32 | ~0.44 | — |

**Panel B — 100 °C, pH 6.8**

| time (min) | ● Glc | ▲ Man | × Fru |
|---|---|---|---|
| 0 | 0 | 0 | 0 |
| 15 | ~0.28 | ~0.20 | ~0.09 |
| 30 | ~0.41 | ~0.28 | ~0.15 |
| 45 | ~0.48 | ~0.31 | ~0.19 |
| 60 | ~0.54 | ~0.33 | ~0.22 |
| 90 | ~0.58 | ~0.35 | ~0.25 |
| 120 | ~0.57 | ~0.36 | ~0.26 |
| 150 | ~0.55 | ~0.37 | ~0.26 |
| 180 | ~0.52 | ~0.36 | ~0.25 |

**Panel C — 120 °C, pH 5.5** (no Fru detected)

| time (min) | ● Glc | ▲ Man | × Fru |
|---|---|---|---|
| 0 | 0 | 0 | not detected |
| 5 | ~0.09 | ~0.13 | — |
| 10 | ~0.21 | ~0.27 | — |
| 15 | ~0.30 | ~0.38 | — |
| 30 | ~0.44 | ~0.56 | — |
| 45 | ~0.49 | ~0.54 | — |
| 60 | ~0.52 | ~0.51 | — |
| 90 | ~0.56 | ~0.47 | — |
| 120 | ~0.57 | ~0.42 | — |

**Panel D — 120 °C, pH 6.8**

| time (min) | ● Glc | ▲ Man | × Fru |
|---|---|---|---|
| 0 | 0 | 0 | 0 |
| 5 | ~0.55 | ~0.30 | ~0.14 |
| 10 | ~0.78 | ~0.42 | ~0.22 |
| 15 | ~0.75 | ~0.40 | ~0.24 |
| 30 | ~0.64 | ~0.35 | ~0.27 |
| 45 | ~0.48 | ~0.31 | ~0.29 |
| 60 | ~0.39 | ~0.27 | ~0.30 |
| 90 | ~0.29 | ~0.22 | ~0.30 |
| 120 | ~0.20 | ~0.18 | ~0.30 |

Authors' statements [M] (§3.1.5, pp. 1656–1658):
> "**Contrary to our expectations, Fru was also observed** in addition to Glc and Man, in samples
> heated at pH initial value 6.8."
> "the formation of Fru showed **no lag phase**, indicating that it can be formed **directly from
> DFG** rather than from isomerisation of Glc and Man."
> "independently of the temperature, at pH 5.5, **Man was formed in higher amounts than Glc**."
> "as the reaction proceeded at 120 °C (Fig. 5C) **Man reached a maximum after 30 min**, decreasing
> afterwards, whereas **Glc did not decrease after 120 min** of heating."
> "At pH 6.8 the opposite was observed, that is, **Glc not only was formed in higher amounts than
> Man but also seemed to be more reactive**, showed by the fast decrease at 120 °C."
> "**Fru, even though formed in considerable amounts was always lower in concentration than the
> other two sugars. Also independently of the reaction conditions no decrease in the Fru amount was
> observed.**"

**Magnitude sanity note, worth carrying downstream:** peak parent-sugar concentrations are
**0.4–0.8 mmol/l on 10 mmol/l DFG, i.e. 4–8 mol%.** Amadori reversal to parent sugars is a
**minority channel** in this system by an order of magnitude. This is the quantitative content
behind the later 2005 Food Chem. 90 conclusion that "the reaction path from DFG into its parents,
glucose and glycine, is not important from a quantitative point of view, even though it does
happen."

### 9.1 The two auxiliary control experiments — RESULTS NOT SHOWN

**Anchor: §3.2.4, p. 1660.** Verbatim:
> "To get a better insight in the sugars enolization pathways two additional experiments were
> performed. **Glc and Fru were heated alone in an aqueous phosphate buffer (0.1 mol/L) solution at
> 120 °C and pH 6.8.** It was observed (**results not shown**) that Glc isomerised preferably into
> Fru, whereas Man was only formed in small amounts. On the other hand, the amount of Glc and Man
> formed from Fru was very low relatively to its decrease. **Fru degraded preferably into acids.**"

And (§3.2.5, p. 1660):
> "These results were compared with the ones obtained from heating Glc and Fru, separately with
> glycine. When the amino acid was present, the **formic acid amount decreased considerably and
> acetic acid became the dominant acid** (**results not shown**). It is worth noting that **no
> glycolic acid or lactic acid were detected**."

**All four of these control experiments are reported with NO numbers whatsoever. Unreadable /
absent by the authors' own choice. They are quotable as qualitative structural evidence only.**

---

## 10. FIGURE 6 — organic acid formation

**Anchor: Fig. 6, p. 1658 (PDF p. 8).** Caption verbatim: "Organic acids formation during thermal
degradation of DFG at 100 (A) and 120 °C (B): acetic acid (○, pH 5.5; ●, pH 6.8); Formic acid
(△, pH 5.5; ▲, pH 6.8)."
Axes as printed: **A** y = `Organic Acids mmol/l` 0–6 (ticks 1), x 0–200; **B** y 0–7 (ticks 1),
x 0–150. Error bars present, definition unstated.

Digitised from PDF p. 8 raster. **All [fig].** Read precision ≈ ±0.15 mmol/l.

**Panel A — 100 °C**

| time (min) | ○ AA pH 5.5 | ● AA pH 6.8 | △ FA pH 5.5 | ▲ FA pH 6.8 |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 15 | ~0.03 | ~1.5 | ~0.02 | ~0.35 |
| 30 | ~0.06 | ~2.5 | ~0.03 | ~0.75 |
| 45 | ~0.10 | ~3.2 | ~0.04 | ~1.05 |
| 60 | ~0.13 | ~3.7 | ~0.05 | ~1.45 |
| 90 | ~0.20 | ~4.4 | ~0.07 | ~1.75 |
| 120 | ~0.28 | ~4.9 | ~0.09 | ~1.95 |
| 150 | ~0.36 | ~4.8 | ~0.12 | ~1.98 |
| 180 | ~0.45 | ~4.8 | ~0.15 | ~2.00 |

**Panel B — 120 °C**

| time (min) | ○ AA pH 5.5 | ● AA pH 6.8 | △ FA pH 5.5 | ▲ FA pH 6.8 |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 5 | ~0.10 | ~1.3 | ~0.06 | ~0.9 |
| 10 | ~0.25 | ~2.9 | ~0.15 | ~1.8 |
| 15 | ~0.40 | ~3.9 | ~0.25 | ~2.3 |
| 30 | ~0.75 | ~4.9 | ~0.50 | ~2.9 |
| 45 | ~1.05 | ~5.2 | ~0.70 | ~3.0 |
| 60 | ~1.25 | ~5.3 | ~0.85 | ~3.05 |
| 90 | ~1.50 | ~5.5 | ~1.00 | ~2.95 |
| 120 | ~1.70 | ~5.6 | ~1.10 | ~2.90 |

Authors' statements [M] (§3.1.6, p. 1658):
> "At pH 5.5 the yield of formic and acetic acid was considerably lower relative to pH 6.8. Also an
> increase of temperature enhanced formation of both acids. **The sum of the yields of both acids
> reached up to 0.8 mmol/mmol of DFG**, which indicates that organic acids, in particular acetic
> acid, are important end products of the Maillard reaction. **Independently of the reaction
> conditions, acetic acid was always formed in higher amounts than formic acid.** The difference
> was more significant at pH 6.8."
> "the level-off observed on organic acids formation (Fig. 6B) matched the complete degradation of
> DFG, which implies **direct involvement of DFG or its early degradation product in the formation
> of organic acids**."

Cross-check [Z]: at 120 °C / pH 6.8, (AA + FA) ≈ 5.6 + 2.9 = 8.5 mmol/l on ~9.8 mmol/l DFG =
**0.87 mmol/mmol**, consistent with the authors' "up to 0.8 mmol/mmol". Acetic acid alone at that
condition ≈ 5.6/9.8 = **57 mol%**, consistent with the abstract's "55 mol%" for acetic acid.

Also [M], on why the pH drop happens (§3.1.6, p. 1659), citing Davidek et al. (ref. 6) [C]:
> "the total amounts of acetic and formic acid were generally lower than the amounts of NaOH added
> to keep the pH constant, suggesting that **other acids might be formed**."

---

## 11. FIGURE 7 — pH drop

**Anchor: Fig. 7, p. 1659 (PDF p. 9).** Caption verbatim: "pH drop during thermal degradation of
DFG at initial reaction pH 6.8. (○, 100 °C; ×, 120 °C)."
Axes as printed: y = `pH`, **6.55–6.85, ticks every 0.05**; x = `Time (min)`, 0–200, ticks 50.
No error bars visible.

Digitised from PDF p. 9 raster. **All [fig].** Read precision ≈ ±0.01 pH unit (the y-axis is
heavily expanded, so this figure digitises unusually well).

| time (min) | ○ 100 °C | × 120 °C |
|---|---|---|
| 0 | ~6.80 | ~6.81 |
| 5 | ~6.79 | ~6.76 |
| 10 | ~6.78 | ~6.72 |
| 15 | ~6.77 | ~6.70 |
| 30 | ~6.74 | ~6.65 |
| 45 | ~6.72 | ~6.62 |
| 60 | ~6.70 | ~6.60 |
| 90 | ~6.67 | ~6.60 |
| 120 | ~6.65 | ~6.60 |
| 150 | ~6.64 | — |
| 180 | ~6.63 | — |

Authors' numbers [M] (§3.1.6, p. 1658):
> "After 120 min at 120 °C the pH drop was more significant than at 100 °C after 180 min, **0.21
> and 0.19**, respectively. Due to the low buffer capacity at initial pH of 5.5 the pH drop for the
> same reaction conditions was **0.8 and 0.9**, respectively, even though the sum of acids formed
> was lower."

Cross-check [Z]: my [fig] read gives 6.81 → 6.60 = **0.21** (120 °C) and 6.80 → 6.63 = **0.17**
(100 °C) vs. the printed 0.19 — within read precision.

**NO FIGURE IS GIVEN FOR THE pH 5.5 DRIFT.** The 0.8 / 0.9 unit drops are text-only, with no time
resolution. **This is the single most important limitation for any pH-axis hold-out built on this
paper: the "pH 5.5" condition is an uncontrolled trajectory from 5.5 down to ~4.6–4.7, and its
shape is unpublished.**

---

## 12. FIGURE 8 — melanoidin formation

**Anchor: Fig. 8, p. 1659 (PDF p. 9).** Caption verbatim: "Melanoidins formation during thermal
degradation of DFG at 100 (△, pH 5.5; ▲, pH 6.8) and 120 °C (○, pH 5.5; ●, pH 6.8)."
Axes as printed: y = `Melanoidins mmol/l`, 0–3, ticks 0.5; x = `Time (min)`, 0–200, ticks 50.
Error bars present, definition unstated.

**Units caution:** by §2.10 this y-axis is **mmol/l of glucose-equivalent sugar incorporated into
the brown polymer**, obtained from A₄₇₀ via ε = 0.64 ± 0.03 l/mmol·cm. It is **not** a melanoidin
molar concentration, and it inherits the ±5% uncertainty of ε.

Digitised from PDF p. 9 raster. **All [fig].** Read precision ≈ ±0.08 mmol/l.

| time (min) | △ 100 °C pH 5.5 | ▲ 100 °C pH 6.8 | ○ 120 °C pH 5.5 | ● 120 °C pH 6.8 |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 5 | ~0.01 | ~0.35 | ~0.10 | ~1.05 |
| 10 | ~0.02 | ~0.75 | ~0.20 | ~1.55 |
| 15 | ~0.03 | ~1.05 | ~0.30 | ~1.85 |
| 30 | ~0.06 | ~1.65 | ~0.55 | ~2.30 |
| 45 | ~0.10 | ~1.95 | ~0.75 | ~2.35 |
| 60 | ~0.15 | ~2.15 | ~0.95 | ~2.40 |
| 90 | ~0.25 | ~2.35 | ~1.30 | ~2.45 |
| 120 | ~0.35 | ~2.45 | ~1.70 | ~2.50 |
| 150 | ~0.48 | ~2.50 | — | — |
| 180 | ~0.60 | ~2.55 | — | — |

Authors' statements [M] (§3.1.7, p. 1659):
> "**pH more than temperature had a strong influence in melanoidins formation from DFG.** Moreover,
> the observed **colour formation is much less than when one starts with Glc/glycine system**, even
> though the amount of DFG formed then is approximately the same as we used as reactant in the
> present study."

That last sentence is the origin of the paper's most consequential structural claim — that DFG is
**not** the dominant browning precursor — which Part II then converts into the finding that 3-DG,
not DFG, is the main colour precursor. **It is quoted here without a supporting number; the
comparison Glc/glycine dataset is cited to ref. 22 (a COST conference proceedings), not shown.**

---

## 13. HETEROCYCLES (HMF, HHMF, DDMP) — a quantitative gap, stated

**Anchor: §3.2.1, pp. 1659–1660.** Verbatim:
> "In the present study **both HHMF and DDMP were identified at both pHs**, which indicates that
> 2,3-enolization occurs not only at pH 6.8 but also at pH 5.5. **Unfortunately, their
> quantification was not possible since no reference material was available.** However, judging by
> the response factor of HMF they were formed in the order of magnitude of **µmoles**. On the other
> hand, **HMF was not formed at pH 6.8 and only µmole amounts at pH 5.5**."

**No numbers. "Order of magnitude of µmoles" is the entire quantitative content.** Against a
10 mmol/l DFG basis this is ≲0.1 mol% — negligible, and the paper says so implicitly by dropping
these species from the Part II model entirely.

---

## NEW-PARAMETER TABLE (consolidated)

**This paper contains NO rate constants, NO equilibrium constants, NO activation energies and NO
fitted parameters.** Every row below is either a measured observable, a measured yield, a
digitised figure value, or a value cited from elsewhere.

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| [DFG]₀ | 10 | mmol (in 100 mL) → 10 mmol/l | all systems | §2.4, p. 1654 | [M] |
| buffer | K₂HPO₄/KH₂PO₄, 0.1 | M | all systems | §2.4, p. 1654 | [M] |
| glycine yield, full grid | see §3 table (43.40–81.46 / 71.02–93.85 / 25.87–74.20 / 83.92–95.82) | % (mol Gly per mol DFG reacted) | 100/120 °C × pH 5.5/6.8, 5–180 min | **Table 1, p. 1657** | **[M]** |
| glycine yield, headline | 83 | mol% | conditions **not stated in the abstract**; best reconstruction pH 6.8 | Abstract, p. 1651 | [M], under-specified |
| acetic acid yield, headline | 55 | mol% | conditions **not stated**; best reconstruction pH 6.8 | Abstract, p. 1651 | [M], under-specified |
| DFG conversion | ~57 (i.e. 57% reacted) | % of initial | 100 °C, pH 5.5, 180 min | §3.1.1, p. 1655 | [M] |
| DFG conversion | 100 ("completely degraded") | % of initial | 120 °C, pH 6.8, 30 min | §3.1.1, p. 1655 | [M] |
| pH/temperature equivalence | −1.3 pH units ≡ +20 °C for DFG degradation rate | — | valid only pH 5–7 | §3.1.1, p. 1655 | [M] |
| pH drop, initial pH 6.8 | 0.19 | pH units | 100 °C, 180 min | §3.1.6, p. 1658 + Fig. 7 | [M] |
| pH drop, initial pH 6.8 | 0.21 | pH units | 120 °C, 120 min | §3.1.6, p. 1658 + Fig. 7 | [M] |
| pH drop, initial pH 5.5 | 0.8 | pH units | 100 °C, 180 min | §3.1.6, p. 1659 | [M] — **no figure, text only** |
| pH drop, initial pH 5.5 | 0.9 | pH units | 120 °C, 120 min | §3.1.6, p. 1659 | [M] — **no figure, text only** |
| DFG time-course | see §5 table | mmol/l | 4 conditions, 0–180 min | Fig. 1, p. 1656 | **[fig]** (PDF p. 6 raster, ±0.3 mmol/l) |
| glycine time-course | see §6 table | mmol/l | 4 conditions, 0–180 min | Fig. 2, p. 1656 | **[fig]** (±0.3 mmol/l) |
| 1-DG and 3-DG time-courses | see §7 tables | mmol/l | 4 conditions, 0–180 min | Fig. 3, p. 1657 | **[fig]** (±0.015 mmol/l — low signal) |
| 3-DG peak | ~0.47 (100 °C) / ~0.43 (120 °C) | mmol/l | pH 5.5 | Fig. 3, p. 1657 | **[fig]** |
| 1-DG peak | ~0.09 (100 °C) / ~0.08 (120 °C) | mmol/l | pH 5.5 | Fig. 3, p. 1657 | **[fig]** |
| 3-DG : 1-DG ratio | 3.2 | dimensionless | 100 °C, pH 6.8 | §3.2.1, p. 1659 | [M] — **contradicted by my [fig] read of ~1.2; see §7 flag** |
| deoxyosones at 5 min | "34%" (denominator not defined) | % | 100 °C, pH 5.5 | §3.1.4, p. 1656 | [M], ambiguous |
| MG time-course | see §8 table | mmol/l | 4 conditions, 0–120 min | Fig. 4, p. 1657 | **[fig]** (±0.08 mmol/l) |
| MG at 120 min | ~1.15 / ~2.00 / ~2.95 / ~3.20 | mmol/l | A/B/C/D as in §8 | Fig. 4, p. 1657 | **[fig]** |
| Glc, Man, Fru time-courses | see §9 tables | mmol/l | 4 conditions | Fig. 5, p. 1658 | **[fig]** (±0.02 mmol/l); **panel↔condition assignment is my reasoned inference, not printed** |
| Glc peak | ~0.78 | mmol/l | 120 °C, pH 6.8, ~10 min | Fig. 5D, p. 1658 | **[fig]** |
| Man peak | ~0.56 | mmol/l | 120 °C, pH 5.5, ~30 min | Fig. 5C, p. 1658 | **[fig]** |
| Fru plateau | ~0.30 | mmol/l | 120 °C, pH 6.8 | Fig. 5D, p. 1658 | **[fig]** |
| Fru | **not detected at all** | — | both temperatures, pH 5.5 | §3.1.5, p. 1657 | [M] — a genuine measured absence |
| total parent-sugar yield | ~4–8 | mol% of [DFG]₀ | all conditions | Fig. 5, p. 1658 | **[Z]** from [fig] |
| AA / FA time-courses | see §10 tables | mmol/l | 4 conditions | Fig. 6, p. 1658 | **[fig]** (±0.15 mmol/l) |
| AA + FA total yield | "up to 0.8" | mmol per mmol DFG | best case (pH 6.8) | §3.1.6, p. 1658 | [M] |
| AA + FA total yield | ~0.87 | mmol/mmol DFG | 120 °C, pH 6.8, 120 min | Fig. 6B, p. 1658 | **[Z]** from [fig] |
| melanoidin time-course | see §12 table | mmol/l (sugar-equivalents) | 4 conditions | Fig. 8, p. 1659 | **[fig]** (±0.08 mmol/l) |
| melanoidin plateau | ~2.5 (pH 6.8) vs ~0.6 (100 °C pH 5.5) / ~1.7 (120 °C pH 5.5) | mmol/l sugar-equiv. | as labelled | Fig. 8, p. 1659 | **[fig]** |
| melanoidin ε at 470 nm | 0.64 ± 0.03 | l/mmol·cm | glucose/glycine | §2.10, p. 1655 (ref. 17) | **[C]** — the primary is `martins2003c.pdf` |
| DFG pKₐ | 8.2 | — | — | §3.1.4, p. 1656 | **[C]** |
| 1-DG/3-DG ratio in a related ARP | 8:5 at pH 4.5; 20:1 at pH 7 | dimensionless | 1-deoxy-1-propylamino-D-fructose, 10 h reflux, phosphate | §3.2.1, p. 1659 (ref. 24, Beck et al. 1988) | **[C]** |
| 3-DG share of α-dicarbonyls | 14 | % of total α-dicarbonyls detected | Glc heated alone, 100 °C, pH 10 | §3.2.5, p. 1660 (ref. 35) | **[C]** |
| Hofmann 3-DG:1-DG | ~4 | dimensionless | Glc/L-alanine, reflux, pH 7 | §3.2.1, p. 1659 (refs. 13, 14) | **[C]** |
| DFG synthesis yield | 19 | % on glycine | in-house prep | §2.2, p. 1653 | [M] |
| DFG melting point | 160–163 | °C | in-house prep | §2.2, p. 1653 | [M] |
| DFG elemental analysis | C 40.13 / H 6.4 / N 5.89 (calcd 40.50 / 6.37 / 5.91) | % | in-house prep | §2.2, p. 1653 | [M] |
| HHMF, DDMP, HMF | "order of magnitude of µmoles"; HMF **absent at pH 6.8** | µmol/l | all conditions | §3.2.1, p. 1660 | [M], **no numeric value published** |
| glycolic acid, lactic acid | **not detected** | — | Glc/Fru + glycine controls | §3.2.5, p. 1660 | [M] |
| Glc-alone and Fru-alone controls | **no numbers published ("results not shown")** | — | 120 °C, pH 6.8, 0.1 mol/L phosphate | §3.2.4, p. 1660 | [M], **unquantified** |
| **rate constants of any kind** | **NONE PRESENT** | — | — | whole paper | negative finding, high confidence |
| **activation energies** | **NONE PRESENT** | — | — | whole paper | negative finding, high confidence |
| **error-bar definition** | **NEVER STATED** | — | all figures | whole paper | defect |

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

**Blocking issue first: Table 1 of this paper is, on my reading, ALREADY DECLARED.** Amendment 1
assigns **FIT** to "Martins 2003 thesis Table 4.1.1 glycine-release yields (65–95%)". That range
and that description match **this paper's Table 1** closely enough that I assess the probability
they are the same measurements at high. **Recommendation: do NOT create a second declaration row
for this paper's Table 1. Instead amend the existing row to name the journal anchor
(`data/articles/martins2003b.pdf`, Table 1, p. 1657) so the two can never be counted twice.** The
existing row's stated range should also be widened to 25.87–95.82% to match the table as printed.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Table 1, columns B + D (pH 6.8)** — glycine yield vs time | **FIT (already declared as thesis Table 4.1.1 — amend, do not duplicate)** | pH | Model-free stoichiometric constraint at the trunk's own pH: it pins the amino-acid release per Amadori consumed without committing to any topology. This is exactly the kind of reactant-side anchor the declaration already reasoned for. |
| **Table 1, columns A + C (pH 5.5)** — glycine yield vs time | **HOLD-OUT** | pH | Genuine off-pH prediction for a pH-fixed trunk, and the *shape* differs qualitatively (non-monotone, dipping to 53.6% at 30 min in column A) not just in level. **Caveat: the pH-5.5 arm drifts 0.8–0.9 pH units mid-run (§11), so a failure here may be a pH-trajectory failure, not a chemistry failure. Score with that stated.** |
| **Fig. 1 (DFG) + Fig. 2 (glycine) + Fig. 6 (acids) + Fig. 8 (melanoidins), pH 6.8 columns only** | **FIT** | pH | These four are the well-resolved, high-signal responses (peak-to-noise ≫ 10) at the trunk's pH. Digitised here at stated precision. |
| **The same four figures, pH 5.5 columns** | **HOLD-OUT** | pH | Orthogonal arm of the same cut. Same drift caveat. |
| **Fig. 5 (parent sugars Glc/Man/Fru), all panels** | **diagnostic_only** | — | Two reasons. (i) The **panel↔condition assignment is my inference, not printed** (§9) — I will not put a fit target on an inferred axis label. (ii) The ION-300 method **co-elutes Man and Fru**, and only the 120-min samples were re-run on the Dionex column that separates them; so most Man points in the underlying dataset are Man + Fru. Structurally this figure is the single most valuable thing in the paper (it is the direct measurement of Amadori → parent-sugar reversal at ~4–8 mol%) — but it should inform the *prior*, not be scored. |
| **Fig. 3 (deoxyosones)** | **diagnostic_only** | — | Signal is 0.005–0.47 mmol/l against a figure resolution of ~0.015; four series overlap; and the authors' own 3-DG:1-DG = 3.2 claim contradicts my digitisation (~1.2). Unsafe as a scored target in either arm. |
| **Fig. 4 (methylglyoxal)** | **FIT for pH 6.8, HOLD-OUT for pH 5.5** | pH | Clean, monotone, well-separated series. Note the authors' "more than doubled" claim overstates the pH 6.8 case (1.6× on my read) — carry that as a flag, not a correction. |
| **Fig. 7 (pH drop at initial pH 6.8)** | **FIT (measured input, not a scored target)** | — | This is a *condition* of the experiment, not an outcome. If the trunk ever models in-run pH, this is the forcing function. It should enter as an input the way the melanoidin ε does, never as something the model is scored on. **The pH 5.5 equivalent does not exist and must be treated as missing, not as zero.** |
| **All figure time-courses collectively** | **DO NOT ALSO ASSIGN VIA PART II** | — | Part II Figs. 2 and 4–7 are the same points. Part II's Table 3 is a regression on them. Assign here once; cross-reference from `martins2003_extraction.md`. Failing to do this is precisely the Cerny-2007-style disjointness violation Amendment 2 had to fix retroactively. |
| **§3.2.4 / §3.2.5 control experiments (Glc alone, Fru alone, ± glycine)** | **not usable — no numbers published** | — | "Results not shown" four times. Quotable as qualitative mechanism only. |
| **HHMF / DDMP / HMF** | **not usable — no numbers published** | — | "quantification was not possible since no reference material was available". |
| **Abstract's "83 and 55 mol%"** | **do not use as a scored target** | — | The abstract does not state the conditions. Use Table 1 (for glycine) and Fig. 6 (for acetic acid) with explicit conditions instead. |

**Cut-axis summary:** the clean cut is **pH (6.8 FIT / 5.5 HOLD-OUT)**, consistent with the
declaration's existing convention for the Martins Ch. 6 pH ladder and with the Hofmann precedent.
Temperature is **not** a viable cut here: with only two levels and no Ea published anywhere in this
paper, holding one temperature out leaves the other unidentifiable.

**Circularity risks, stated plainly:**
1. **Highest risk — the Table 1 / thesis Table 4.1.1 identity.** If a prior wave already fit thesis
   Table 4.1.1 and a new wave now fits this paper's Table 1 as if new, the same measurements enter
   the fit twice and the apparent evidence base inflates with no new information.
2. **Second risk — Part I / Part II.** Any hold-out built from these figures is contaminated if
   Part II's Table 3 rate constants are fit, because Table 3 *was regressed on these exact points*.
   The two cannot both be independent. **Recommendation: pick one side. Either fit Part II's Table 3
   (parameter-space anchor) and hold out the pH 5.5 time-courses, or fit the pH 6.8 time-courses
   directly and treat Table 3 as diagnostic. Not both.**
3. **Third risk — melanoidins.** The melanoidin y-axis in Fig. 8 is A₄₇₀ divided by the ε from
   `martins2003c.pdf`. The declaration already assigns that ε as "FIT (measured input, not a scored
   target)". Consistent — but it means Fig. 8 is **not** an independent browning measurement; it is
   an absorbance measurement wearing a concentration's units, and its uncertainty is the
   convolution of A₄₇₀ scatter with ε's ±5%.
