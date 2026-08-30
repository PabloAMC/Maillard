# Bagiyan, Koroleva, Soroka & Ufimtsev 2004 (Kinet. Catal. 45(3) 372–380; no DOI printed) — Wave K4a extraction 2026-08-28

**Source PDF:** `data/articles/bagiyan2004.pdf` (9 pp.). Read method: both — born-digital text layer
(`pdftotext -layout`) plus 200-dpi page rasters of PDF pp. 3, 5, 6, 7 to verify Tables 1–4 cell by
cell (the text layer scrambled the Table 1 pH-10/11/13 blocks and the Table 4 "Notes" column, both
now corrected from the rasters).

## 0. PAPER IDENTITY — PARTIALLY MATCHES; **THE EXPECTED STRUCTURE–REACTIVITY TABLE WITH THIOL pKa VALUES IS NOT IN THIS PAPER**

| field | value |
|---|---|
| authors | G. A. Bagiyan, I. K. Koroleva, N. V. Soroka, A. V. Ufimtsev |
| title | "Kinetics of the Catalytic Oxidation Reactions of Thiol Compounds in Aqueous Solutions in the Presence of Copper Ions" |
| venue | *Kinetics and Catalysis*, Vol. 45, No. 3, 2004, pp. 372–380. Translated from *Kinetika i Kataliz*, Vol. 45, No. 3, 2004, pp. 398–406. Copyright © 2004 MAIK "Nauka/Interperiodica". Article code 0023-1584/04/4503-0372 |
| DOI | **not printed in the PDF** |
| received | Received August 20, 2002 |
| affiliation | Konstantinov Institute of Nuclear Physics, Russian Academy of Sciences, Gatchina, Leningrad oblast, 188300 Russia |
| PDF character | born-digital (FrameMaker5 → Acrobat Distiller 4.05), clean text layer; English translation of a Russian original |

**The file is the right paper by title and subject — Cu(I)-catalysed O₂ oxidation of thiols, kinetic
tables, pH-resolved — but it does NOT contain the artefacts the wave brief expected.** Specifically:

- **There is no table of thiol pKa values.** The word pKa does not appear in the paper. The only
  reference to thiol dissociation is qualitative: "an ascending branch of this pH dependence most
  frequently but not always replicates a dissociation curve of the –SH groups of thiol compounds"
  (p. 376).
- **There are no rate constants k anywhere in the paper.** Every tabulated kinetic number is an
  **initial rate w₀ in mol l⁻¹ s⁻¹**, not a rate constant. The rate law (Eq. 3) has variable and
  fractional orders and the authors explicitly refuse to treat it quantitatively (see §1).
- **There is no relative-rate column and no fitted structure–reactivity correlation.** The only
  structure–reactivity content is two *qualitative* reactivity orderings written out in prose (§5)
  and one prose ratio (250-fold, homocysteine vs butylmercaptan) with no supporting table.
- The **only** regression-like relation the authors fit is a Brønsted–Bjerrum ionic-strength
  relation (Eq. 6) for one compound.

If the wave expected the Bagiyan structure–reactivity dataset **with pKa values**, that is a
different paper (this one's own refs. 6 and 11, Bagiyan et al., *Izv. Akad. Nauk, Ser. Khim.* 2002,
no. 5, pp. 1069 and 1075, are the companion stoichiometry/kinetics papers). **Flagged for the
orchestrator: a wrong-file / wrong-expectation call, not a wrong download — the title matches.**

**Nomenclature defect to carry forward.** The table titles and table bodies write the compounds as
**TSN, ESN, DSN**, while the running text writes the same compounds as **TSH, ESH, DSH**. This is a
Cyrillic-Н → Latin-N transliteration artefact surviving the translation. TSN ≡ TSH ≡
(CH₃)₃N⁺CH₂CH₂SH; ESN ≡ ESH ≡ NH₂CH₂CH₂SH; DSN ≡ DSH ≡ (CH₃)₂NCH₂CH₂SH. Transcribed below exactly as
printed in each location.

## 1. ONE-PARAGRAPH VERDICT

This paper gives the model **no rate constants, no pKa values, no activation energies and no
temperature series** — every kinetic experiment is at a single temperature, 20 °C. What it does give
is a large, dense, pH-resolved matrix of **initial O₂-consumption rates w₀ (mol l⁻¹ s⁻¹)** for three
representative thiols (TSH, ESH, DSH) across pH 6–13 at varying [thiol]₀ and [Cu⁺], plus **reaction
orders x (in [O₂]), y (in [RSH]) and z (in [Cu⁺])** at each pH, and a fourth table on the ESH-induced
acceleration of TSH oxidation. The authors state explicitly that these orders are variable and
fractional and that the rate law "can only serve as the qualitative characteristics of catalytic
oxidation reactions" — i.e. **the authors themselves disclaim a quantitative rate constant.** The
paper's usable content for a Maillard/thiol model is therefore *mechanistic and ordinal*: (a) the
kinetic classification of thiols into three groups by chelating ability; (b) two prose reactivity
orderings spanning ~13 thiols; (c) the pH range over which Cu-catalysed thiol loss switches from
O₂-first-order to O₂-zero-order; (d) the observation that ionic strength has no effect on
nonchelating thiols but accelerates chelating aminothiols in alkali; (e) the 250-fold
homocysteine-vs-butylmercaptan spread at equal [Cu⁺]. There is also one Brønsted-type slope
(Eq. 6, coefficient 1.02) and a set of τ½ values for O₂ binding (Fig. 6).

## 2. SYSTEM DEFINITION — verbatim

| condition variable | value as printed |
|---|---|
| temperature | **20 °C** (stated in the title line of Tables 1, 2 and 3; the only temperature in the paper). Cuvettes are "closed thermostatted 3- to 10-ml cuvettes" |
| buffer | **phosphate–borate buffer solutions (0.1–0.3 mol/l)**; "initial components of high-purity grade were used". For the ionic-strength series the buffer contained "sodium orthophosphate and sodium tetraborate ((6.25–62.5) × 10⁻³ mol/l each)" |
| pH | studied "in solutions with pH 5–13"; tabulated at pH 6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0 (Table 2), 6.0–11.0 and 13.0 (Table 1), 6.0–11.0 (Table 3), 11 (Table 4), 10.5 (Figs. 4, 5), 8.5/11.0/11.5 (Fig. 6), 7.5/10.0 (Fig. 1) |
| hot-pH correction | not applicable (20 °C); not discussed |
| [aminothiol]₀ | "varied within the range 1 × 10⁻⁴ – 1 × 10⁻¹ mol/l" |
| [Cu⁺] | "the concentration of copper as the acetonitrile complex **[Cu(CH₃CN)₄]ClO₄** was varied within the range 1 × 10⁻⁷ – 1 × 10⁻⁴ mol/l"; expressed in the tables as **g-ion/l** |
| [O₂]₀ | **3 × 10⁻⁴ mol/l** in every table and figure caption |
| atmosphere / headspace | closed cell; "To prevent atmospheric oxygen from entering the system, the cell was placed in closed thermostatted 3- to 10-ml cuvettes with test solutions with the use of ground-glass joints"; reagents supplied as "deaerated concentrated solutions" through polyethylene capillaries |
| volume | 3–10 ml cuvettes |
| agitation | "A magnetic stirrer with a constant speed of rotation was used for stirring" |
| quench | none — continuous polarographic O₂ monitoring, "in the automatic mode from changes in the cathode current intensity (a signal from the polarographic cell was applied to the input of an EPP-09 recorder) rather than with the use of a sampling method" |
| analytical method | modified membrane-covered polarographic O₂ electrode: Pt cathode vs Ag/AgCl reference, potential difference **0.6 V**; "Providing a gap of about 2 µm between the cathode and the membrane in the design of a membrane electrode and using Teflon films 10–15 µm in thickness as membranes, we shortened the time taken to reach a steady state at the cathode to **2–3 s**" (vs 40–60 s for a standard Clark electrode) |
| measurement error | "This allowed us to detect the consumption of O₂ in solution for several minutes with an error of **5% or lower**"; "The measurement error in the concentration of O₂ in solutions was no higher than **5%** at [O₂] = (2–100) × 10⁻⁵ mol/l" |
| replication (n) | no per-point n. Aggregate only: "this procedure allowed us to determine the kinetic characteristics of **30 thiol compounds** for a great body of kinetic experiments (**more than 800**)"; abstract says "more than 30 thiols" |
| error bars | **none printed anywhere.** No ± values, no SE, no R² in any table |
| thiol purity | "Aminothiols were synthesized and purified to remove possible impurities in accordance with a published procedure [9]" (Bagiyan et al., *Zh. Obshch. Khim.* 1976, 46, 365) |
| jet-mixing experiments (Fig. 6) | phosphate–borate buffer with O₂ 3 × 10⁻⁴ mol/l mixed with aminothiol + Cu⁺ at a **20 : 1 flow ratio**, giving resulting [RSH] = (1–0.5) × 10⁻² mol/l and [Cu⁺] = (1–2) × 10⁻⁴ mol/l; **dead time ~30 ms** |

**Order determination, verbatim (p. 373–374):** initial-rate ratios,
`(d[O₂]'/dt')/(d[O₂]''/dt'') = ([O₂]'/[O₂]'')ⁿ` (Eq. 1) and
`n = {log(d[O₂]'/dt') − log(d[O₂]''/dt'')}/(log[O₂]' − log[O₂]'')` (Eq. 2),
"where n is the order of reaction with respect to O₂. The orders of reaction with respect to thiol
compounds and Cu⁺ were calculated from [O₂] = f(t) functions found at various concentrations of thiol
compounds and Cu⁺ with the use of equations analogous to Eqs. (1) and (2)."

**The authors' own disclaimer, verbatim (p. 374), which governs how any of these numbers may be used:**
> "Note that, under conditions of variable and fractional orders, rate equations of the type
> −d[O₂]/dt = w₀ = k[O₂]ˣ[RSH]ʸ[Cu⁺]ᶻ  (3)
> can only serve as the qualitative characteristics of catalytic oxidation reactions. Only rate
> equations in which the orders of reaction with respect to reactants are constant over the tested
> range of concentrations (pH) can be correct for quantitative kinetic analysis."

The symbol **k in Eq. (3) is never evaluated anywhere in the paper.**

## 3. TABLE 1 — TSN (trimethylammonium-ethanethiol), group 1 (nonchelating / weakly chelating)

**Anchor: Table 1, p. 374 (PDF page 3).** Title as printed:
"Table 1. Kinetic data on the catalytic reaction of O₂ with (CH₃)₃N⁺CH₂CH₂SH (TSN) in the presence of
copper ions (phosphate–borate buffer solutions, 20°C). [O₂]₀ = 3 × 10⁻⁴ mol/l"

Column headers as printed (the table is set as two side-by-side half-tables sharing one header):
`pH` | `[TSN]₀ × 10³, mol/l` | `[Cu⁺] × 10⁶, g-ion/l` | `w₀ × 10⁷, mol l⁻¹ s⁻¹` | `Order of reaction*`
Footnote as printed: `* See Eq. (3).`

Raster-verified against PDF p. 3. Transcribed in printed reading order, left half-table then right
half-table. All values [M].

**Left half-table**

| pH | [TSN]₀ × 10³, mol/l | [Cu⁺] × 10⁶, g-ion/l | w₀ × 10⁷, mol l⁻¹ s⁻¹ | Order of reaction* |
|---|---|---|---|---|
| 6.0 | 2.2 | 5.0 | 0.18 | x = 1.0 |
|  | 2.2 | 50.0 | 0.30 | z = 0.22 |
|  | 100 | 5.0 | 1.35 | y = 0.73 |
|  | 100 | 50.0 | 4.85 | z = 0.55 |
| 8.0 | 2.2 | 5.0 | 0.82 | x = 1.0 |
|  | 2.2 | 50.0 | 4.15 | z = 0.59 |
|  | 20.0 | 25.0 | 2.65 | y = 0.25 |
|  | 20.0 | 50.0 | 5.35 | x = 1.0 |
| 10.0 | 2.2 | 5.0 | 0.70 | x = 1.00 |
| 11.0 | 2.2 | 5.0 | 0.70 | x = 1.00 |
|  | 2.2 | 50.0 | 6.50 | z = 0.95 |
|  | 20.0 | 25.0 | 1.50 | y = 0.11 |
|  | 20.0 | 50.0 | 3.00 | z = 1.00 |

**Right half-table**

| pH | [TSN]₀ × 10³, mol/l | [Cu⁺] × 10⁶, g-ion/l | w₀ × 10⁷, mol l⁻¹ s⁻¹ | Order of reaction* |
|---|---|---|---|---|
| 7.0 | 2.2 | 5.0 | 0.54 | x = 1.0 |
|  | 2.2 | 50.0 | 0.95 | z = 0.25 |
|  | 110 | 25.0 | 2.70 | y = 0.50 |
|  | 110 | 50.0 | 4.50 | z = 0.75 |
| 9.0 | 2.2 | 5.0 | 0.62 | x = 1.0 |
|  | 2.2 | 50.0 | 4.15 | z = 0.83 |
|  | 20.0 | 25.0 | 2.10 | y = 0.03 |
|  | 20.0 | 50.0 | 4.10 | z = 1.00 |
| 10.0 | 2.2 | 50.0 | 6.50 | z = 0.95 |
| 13.0 | 2.2 | 50.0 | 3.35 | x = 1.00 |
|  | 2.2 | 50.0 | 1.68 | y = 0.63 |
|  | 20.0 | 25.0 | 0.84 | z = 1.00 |

*Printing anomaly, reported not fixed:* the pH-13.0 block lists two rows with **identical**
[TSN]₀ = 2.2 and [Cu⁺] = 50.0 but different w₀ (3.35 and 1.68). As printed this is internally
inconsistent; raster-checked twice, the cells read as transcribed. Also the pH-8.0 left block prints
`x = 1.0` twice within the same four-row group (rows 1 and 4), which cannot both be a fresh
determination for that group.

## 4. TABLE 2 — ESN (2-aminoethanethiol / cysteamine), group 2 (chelating)

**Anchor: Table 2, p. 376 (PDF page 5).** Title as printed:
"Table 2. Kinetic data on the catalytic reaction of O₂ with NH₂CH₂CH₂SH (ESN) in the presence of
copper ions (phosphate–borate buffer solutions, 20°C). [O₂]₀ = 3 × 10⁻4 mol/l"

Column headers as printed: `pH` | `[ESN]₀ × 10³, mol/l` | `[Cu⁺] × 10⁶, g-ion/l` | `w₀ × 10⁷, mol l⁻¹ s⁻¹` | `Order of reaction*`. Footnote: `* See Eq. (3).`
(Note the title's exponent is set as `10-4` without the superscript minus; read as 3 × 10⁻⁴ mol/l.)

Raster-verified against PDF p. 5. All values [M].

**Left half-table**

| pH | [ESN]₀ × 10³, mol/l | [Cu⁺] × 10⁶, g-ion/l | w₀ × 10⁷, mol l⁻¹ s⁻¹ | Order of reaction* |
|---|---|---|---|---|
| 6.0 | 15.0 | 500.0 | 8.10 | x = 1.00 |
|  | 100.0 | 5.0 | 1.35 | y = 0.60 |
|  | 100.0 | 50.0 | 2.40 | z = 0.25 |
|  | 2.3 | 5.0 | 0.15 | – |
|  | 2.3 | 50.0 | 0.24 | – |
| 7.0 | 10.0 | 1.0 | 0.99 | – |
| 7.5 | 100.0 | 5.0 | 12.9 | x = 1.00 |
|  | 100.0 | 25.0 | 34.0 | z = 0.52 |
| 8.0 | 100.0 | 5.0 | 9.10 | x = 0.90 |
|  | 100.0 | 25.0 | 49.0 | z = 1.00 |
|  | 2.2 | 5.0 | 7.50 | y = 0.42 |
|  | 2.2 | 25.0 | 19.0 | z = 0.60 |
| 10.0 | 100.0 | 5.0 | 12.0 | x = 0.35 |
|  | 100.0 | 25.0 | 170.0 | z = 1.75 |
|  | 2.2 | 5.0 | 46.0 | y = 0.35 |
|  | 2.2 | 25.0 | 510 | z = 1.50 |
| 12.0 | 100.0 | 5.0 | 1.55 | x = 0 |
|  | 100.0 | 25.0 | 42.0 | z = 2.00 |
|  | 2.2 | 5.0 | 1.60 | y = 0 |
|  | 2.2 | 25.0 | 43.0 | z = 2.00 |

**Right half-table**

| pH | [ESN]₀ × 10³, mol/l | [Cu⁺] × 10⁶, g-ion/l | w₀ × 10⁷, mol l⁻¹ s⁻¹ | Order of reaction* |
|---|---|---|---|---|
| 7.0 | 10.0 | 50.0 | 2.70 | x = 1.00 |
|  | 100.0 | 50.0 | 36.0 | z = 0.30 |
|  | 1.0 | 50.0 | 0.60 | y = 0.81 |
|  | 100.0 | 5.0 | 12.0 | – |
|  | 2.2 | 50.0 | 1.11 | – |
|  | 2.2 | 5.0 | 0.54 | – |
| 7.5 | 2.2 | 5.0 | 2.10 | y = 0.48 |
|  | 2.2 | 50.0 | 6.50 | – |
| 9.0 | 100.0 | 25.0 | 69.0 | x = 0.75 |
|  | 100.0 | 5.0 | 10.5 | z = 1.25 |
|  | 5.0 | 5.0 | 41.0 | y = 0.36 |
|  | 2.2 | 2.5 | 20.7 | z = 1.00 |
| 11.0 | 100.0 | 25.0 | 16.2 | x = 0.10 |
|  | 100.0 | 5.0 | 0.65 | z = 2.00 |
|  | 2.2 | 25.0 | 596 | y = 0.34 |
|  | 2.2 | 5.0 | 2.85 | z = 1.90 |
| 13.0 | 100.0 | 25.0 | 45.0 | x = 0 |
|  | 100.0 | 5.0 | 3.20 | z = 1.50 |
|  | 2.2 | 25.0 | 47.0 | y = 0 |
|  | 2.2 | 5.0 | 3.20 | z = 1.50 |

*Anomaly, reported not fixed:* at pH 11.0, [ESN]₀ = 100 mmol/l, w₀ falls from 16.2 (25 µM Cu⁺) to
0.65 (5 µM Cu⁺) — a 25-fold drop for a 5-fold [Cu⁺] change, consistent with the stated z = 2.00; but
at [ESN]₀ = 2.2 mmol/l the same pH gives 596 vs 2.85, a 209-fold drop for the same 5-fold [Cu⁺]
change, which is z ≈ 3.3, not the stated z = 1.90. Raster-verified as printed.

## 5. TABLE 3 — DSN (dimethylaminoethanethiol), group 3 (intermediate / weakly chelating)

**Anchor: Table 3, p. 377 (PDF page 6).** Title as printed:
"Table 3. Kinetic data on the catalytic reaction of O₂ with (CH₃)₂NCH₂CH₂SH (DSN) in the presence of
copper ions (phosphate–borate buffer solutions, 20°C). [O₂]₀ = 3 × 10⁻⁴ mol/l"

Column headers as printed: `pH` | `[DSN]₀ × 10³, mol/l` | `[Cu⁺] × 10⁶, g-ion/l` | `w₀ × 10⁷, mol l⁻¹ s⁻¹` | `Order of reaction*`. Footnote: `* See Eq. (3).`

Raster-verified against PDF p. 6. All values [M].

**Left half-table**

| pH | [DSN]₀ × 10³, mol/l | [Cu⁺] × 10⁶, g-ion/l | w₀ × 10⁷, mol l⁻¹ s⁻¹ | Order of reaction* |
|---|---|---|---|---|
| 6.0 | 10.0 | 6.6 | 1.30 | x = 1.00 |
|  | 10.0 | 84.5 | 1.70 | z = 0.10 |
|  | 5.0 | 6.6 | 0.60 | y = 1.00 |
| 8.0 | 10.0 | 6.6 | 16.5 | x = 1.00 |
|  | 10.0 | 84.5 | 52.0 | z = 0.45 |
|  | 5.0 | 84.5 | 32.0 | y = 0.68 |
| 10.0 | 10.0 | 6.6 | 33.0 | x = 0.70 |
|  | 10.0 | 84.5 | 430 | z = 1.00 |
|  | 5.0 | 6.6 | 31.0 | y = 0 |

**Right half-table**

| pH | [DSN]₀ × 10³, mol/l | [Cu⁺] × 10⁶, g-ion/l | w₀ × 10⁷, mol l⁻¹ s⁻¹ | Order of reaction* |
|---|---|---|---|---|
| 7.0 | 5.0 | 6.6 | 3.00 | x = 1.00 |
|  | 10.0 | 6.6 | 6.30 | z = 1.00 |
|  | 10.0 | 84.5 | 9.00 | y = 0.14 |
| 9.0 | 10.0 | 6.6 | 24.0 | x = 1.00 |
|  | 10.0 | 84.5 | 220 | z = 0.87 |
|  | 5.0 | 84.5 | 175 | y = 0.30 |
| 11.0 | 10.0 | 6.6 | 46.0 | x = 0 |
|  | 10.0 | 25.0 | 270 | z = 1.30 |
|  | 5.0 | 6.6 | 50.0 | y = 0.12 |

## 6. TABLE 4 — mixed-thiol system: TSH accelerated by ESH at pH 11

**Anchor: Table 4, p. 378 (PDF page 7).** Title as printed:
"Table 4. Kinetics of the catalytic oxidation reaction of TSH in an alkaline medium (pH 11) in the
presence of ESH. [O₂]₀ = 3 × 10⁻⁴ mol/l"

Column headers as printed: `[TSN]₀ × 10³, mol/l` | `[ESN]₀ × 10³, mol/l` | `[Cu⁺] × 10⁵, g-ion/l` |
`w₀ × 10⁷, mol l⁻¹ s⁻¹` | `Notes`
(Note: the title says TSH/ESH, the column heads say TSN/ESN — the transliteration defect again.
Note also that [Cu⁺] here is **× 10⁵**, not × 10⁶ as in Tables 1–3.)

Raster-verified against PDF p. 7. **The "Notes" entries are block-level, spanning the horizontal rule
groups shown below — they are NOT one-per-row.** The text layer misaligned them; this is the raster
reading. All values [M] except the Notes, which are [F] (fitted orders).

**Block A (10 rows) — Notes for the block: `x = 1.00`, `z = 1.00`, `y_ESN = 0.80`, `y_TSN = 0.63`**

| [TSN]₀ × 10³, mol/l | [ESN]₀ × 10³, mol/l | [Cu⁺] × 10⁵, g-ion/l | w₀ × 10⁷, mol l⁻¹ s⁻¹ |
|---|---|---|---|
| 10.0 | – | 5.0 | 3.4 |
| 10.0 | 0.2 | 5.0 | 0.4 |
| 10.0 | – | 2.5 | 1.4 |
| 10.0 | 0.2 | 2.5 | 4.5 |
| 10.0 | 1.0 | 2.5 | 17.0 |
| 10.0 | 1.0 | 5.0 | 35.0 |
| 50.0 | 0.2 | 5.0 | 4.1 |
| 50.0 | 1.0 | 5.0 | 9.4 |
| 50.0 | 1.0 | 10.0 | 17.5 |
| 50.0 | – | 5.0 | 3.3 |

**Block B (8 rows) — Notes for the block, in printed order: `x = 0`, `z = 2.00`, `y_TSN = 1.00`,
`z = 1.90`, `y_ESN = 0`, `z = 1.25`, `x = 0.65`, `x = 1.00`, `z = 1.10`**

| [TSN]₀ × 10³, mol/l | [ESN]₀ × 10³, mol/l | [Cu⁺] × 10⁵, g-ion/l | w₀ × 10⁷, mol l⁻¹ s⁻¹ |
|---|---|---|---|
| 5.0 | 5.0 | 1.0 | 12.0 |
| 5.0 | 5.0 | 0.5 | 1.4 |
| 10.0 | 5.0 | 0.5 | 2.7 |
| 10.0 | 5.0 | 1.0 | 12.0 |
| 10.0 | 2.5 | 1.0 | 12.0 |
| 10.0 | 2.5 | 0.5 | 4.8 |
| 10.0 | 1.8 | 0.5 | 8.4 |
| 10.0 | 1.8 | 1.0 | 19.5 |

*Anomaly, reported not fixed:* Block B has **8 data rows but 9 Notes entries**, so the mapping of
orders to rows is not recoverable from the print. Also, block A row 2 shows ESH *suppressing* TSH
oxidation (3.4 → 0.4 at 5.0 × 10⁻⁵ Cu⁺) while row 4 shows the same 0.2 mmol/l ESH *accelerating* it
(1.4 → 4.5 at 2.5 × 10⁻⁵ Cu⁺) — the paper's narrative is entirely about acceleration and never
addresses the 3.4 → 0.4 row.

## 7. STRUCTURE–REACTIVITY CONTENT — prose orderings only, no table, no fit

**Anchor: p. 374 (PDF page 3), group 1.** Verbatim:
> "Kinetic behavior of this kind was found for the following thiol compounds, which are arranged in
> order of decreasing Cu⁺ activity in the reactions of their oxidation:
> NH₂CH(COOH)CH₂CH₂SH (HcySH) > NH₂(CH₂)₃SH (PSH) > morpholinoethanethiol (MfSH) >
> piperidinoethanethiol (PipSH) > NH₂(CH₂)₄SH (BSH) > NH₂(CH₂)₅SH (PentSH) >
> (CH₃)₃N⁺CH₂CH₂SH (TSH) > C₂H₅SH (EtSH) > C₃H₇SH (PrSH) > C₄H₉SH (BuSH)."

Rate law for this group (Eq. 4, p. 374): `w₀ = k[O₂]¹[RSH]^(1→0)[Cu⁺]^(0→1)` — orders written as
pH-dependent ranges, not numbers.

**Anchor: p. 374 (PDF page 3), group 2.** Verbatim:
> "This oxidation behavior was typical of NH₂CH₂CH₂SH (ESH) (Table 2) and the following thiol
> compounds:
> NH₂CH(C₂H₅)CH₂SH ≈ NH₂CH(CH₃)CH₂SH > (CH₃)NHCH₂CH₂SH > NH₂CH₂CH₂SH (ESH) >
> NH₂CH(COOH)CH₂SH (**CySH**) > HS-CH₂-COOH (TGA) > glutathione (**GSH**) > NH₂C(C₂H₅)CH₂SH >
> NH₂CH₂CH₂CH₂NHCH₂CH₂SH."

Rate law (Eq. 5, p. 374): `w₀ = k[O₂]^(1→0)[RSH]^(1→0)[Cu⁺]^(0→2)`.

**Anchor: p. 375 (PDF page 4), group 3.** Verbatim:
> "The catalytic oxidation reactions of (CH₃)₂NCH₂CH₂SH (DSH) (Table 3) and the following thiol
> compounds occur in this manner:
> (C₂H₅)₂NCH₂CH₂SH > (CH₃)₂NCH₂CH₂SH (DSH) > NH₂CH₂C(CH₃)₂SH > HOCH₂CH₂CH₂SH > HOCH₂CH₂SH.
> Weakly chelating thiol compounds belong to this group."

**The only quantitative structure–reactivity statement in the whole paper**, anchor p. 375 (PDF
page 4), verbatim:
> "Thus, the rate of O₂ consumption in the catalytic oxidation reaction of homocysteine was higher
> than that in the catalytic oxidation reaction of butylmercaptan by a factor of **250** at equal
> concentrations of Cu⁺ in the solutions."

**No table anywhere in the paper supports this 250-fold figure** — neither HcySH nor BuSH appears in
Tables 1–4. Provenance [M] as reported, but with no traceable underlying data in this file. The pH,
[thiol]₀ and [Cu⁺] at which the factor of 250 holds are not stated.

**Cysteine and glutathione — the two thiols a Maillard model would most want — appear only as names
in the group-2 ordering.** No w₀, no rate constant and no pKa is given for either. CySH appears
additionally as the subject of Fig. 3 (see §8).

## 8. FIGURES — captions transcribed; conditions only, no digitisable rate table

| figure | anchor | caption content (verbatim conditions) | numbers extractable |
|---|---|---|---|
| Fig. 1 | p. 373 (PDF p. 2) | "The time dependence of [O₂] in the catalytic oxidation of **cysteine** in the presence of Cu⁺ ions: (a) pH 7.5; (1) [thiol]₀ = 2.5 × 10⁻³ mol/l; [Cu⁺] = 4.0 × 10⁻⁷ g-ion/l; (2) [thiol]₀ = 5 × 10⁻³ mol/l; [Cu⁺] = 4.0 × 10⁻⁷ g-ion/l; (3) [thiol]₀ = 5 × 10⁻³ mol/l; [Cu⁺] = 8.0 × 10⁻⁷ g-ion/l; [O₂]₀ = 3 × 10⁻⁴ mol/l; (b) pH 10.0; (1) [thiol]₀ = 2.5 × 10⁻³ mol/l; [Cu⁺] = 1.2 × 10⁻⁶ g-ion/l; (2) [thiol]₀ = 5 × 10⁻³ mol/l; [Cu⁺] = 1.2 × 10⁻⁶ g-ion/l; (3) [thiol]₀ = 2.5 × 10⁻³ mol/l; [Cu⁺] = 2.4 × 10⁻⁶ g-ion/l; [O₂]₀ = 3 × 10⁻⁴ mol/l." | **This is the ONLY cysteine kinetic data in the paper, and it is figure-only [fig].** Axes: [O₂] × 10⁴ mol/l, 0–2.0+; τ, 0–20 min. I have not digitised it. |
| Fig. 2 | p. 375 (PDF p. 4) | "The pH dependence of the rates of oxidation of nonchelating and weakly chelating thiol compounds by molecular oxygen at [thiol]₀ = 5 × 10⁻³ mol/l and [O₂]₀ = 3 × 10⁻⁴ mol/l. (1, 2) NH₂CH(COOH)(CH₂)₂SH and (3) (CH₃)₃NCH₂CH₂SH. [Cu⁺], g-ion/l: (1) 2.5 × 10⁻⁶, (2) 2.5 × 10⁻⁵, and (3) 5.0 × 10⁻⁵." | log w₀ axis −7 to −5; pH axis 8–12. **Homocysteine (HcySH) pH profile, figure-only [fig].** Not digitised. |
| Fig. 3 | p. 375 (PDF p. 4) | "The pH dependence of the rates of oxidation of the chelating thiol compound NH₂CH(COOH)CH₂SH by molecular oxygen at [thiol]₀ = 5 × 10⁻³ mol/l and [O₂]₀ = 3 × 10⁻⁴ mol/l. [Cu⁺], g-ion/l: (1) 5.0 × 10⁻⁶ and (2) 5.0 × 10⁻⁵." | **This is CYSTEINE (CySH).** log w₀ axis −7 to −5; pH axis 7–11; bell-shaped. Figure-only [fig], not digitised. |
| Fig. 4 | p. 377 (PDF p. 6) | "Dependence of the rate of (CH₃)HNCH₂CH₂SH oxidation on the ionic strength of a phosphate–borate buffer solution. [Thiol]₀ = 1 × 10⁻³ mol/l; [Cu⁺] = 5.0 × 10⁻⁶ g-ion/l. The ionic strength was adjusted with phosphate–borate buffer solutions containing sodium orthophosphate and sodium tetraborate ((6.25–62.5) × 10⁻³ mol/l each) [O₂]₀ = 3 × 10⁻⁴ mol/l." **Note: the compound in this caption, (CH₃)HNCH₂CH₂SH, is N-methylcysteamine — NOT the DSN of Table 3 on the same page.** | log w₀ axis −7.0 to −6.2; x-axis √I from ~0.1 to ~0.55. Linear, ~6 points, [fig]. Not digitised. |
| Fig. 5 | p. 378 (PDF p. 7) | "The time dependence of [O₂] in the catalytic oxidation of ethylmercaptan: (1) without the addition of ESSE and (2) [ESSE] = 2 × 10⁻⁴ mol/l. [CuSO₄] = 5 × 10⁻⁵ mol/l; [C₂H₅SH]₀ = 1 × 10⁻² mol/l; pH 10.5; [O₂]₀ = 3 × 10⁻⁴ mol/l." **Note: this is the only place CuSO₄ (i.e. Cu²⁺) is used rather than the Cu(I) acetonitrile complex.** | [fig], not digitised |
| Fig. 6 | p. 379 (PDF p. 8) | "Binding of O₂ with the complexes of aminothiols + Cu⁺ for chelating and nonchelating aminothiols: (1) [ESH] = 5 × 10⁻² mol/l; [Cu⁺] = 1 × 10⁻⁴ g-ion/l; pH 8.5; **τ½ = 2.90 s**; (2) [ESH] = 5 × 10⁻² mol/l; [Cu⁺] = 1 × 10⁻⁴ g-ion/l; pH 11.5; **τ½ = 0.78 s**; (3) [PSH] = 5 × 10⁻² mol/l; [Cu⁺] = 2 × 10⁻⁴ g-ion/l; pH 11.0; **τ½ = 0.29 s**." | **Three τ½ values printed in the caption — these are the only half-lives in the paper.** [M]. |

## 9. THE ONE FITTED RELATION IN THE PAPER

**Anchor: Eq. (6), p. 376 (PDF page 5).** As printed:
`log w'₀ = log w₀ + 1.02 Z_a Z_b √I`   (6)
"where I = 1/2 ΣZ²C for the interaction of particles with Z_a = Z_b = −1 [14], where [RS⁻]/[Cu⁺] =
100 : 1 under experimental conditions (pH 10.5) (Fig. 4)."

| parameter | value | provenance |
|---|---|---|
| Brønsted–Bjerrum coefficient | 1.02 | [F] (the only fitted coefficient in the paper); note ref. [14] is Laidler, *Reaction Kinetics*, so the 1.02 is a textbook constant applied, arguably [C] |
| Z_a Z_b | (−1)(−1) = +1 | [M] (stated) |
| ionic-strength range with **no effect** on nonchelating thiols | µ = 0.1–0.5 mol/l | [M] |
| ionic-strength range where chelating aminothiols **accelerate** (alkaline medium) | µ = 0.1–0.5 mol/l | [M] |

**Non-effects, verbatim (p. 376), all [M], all useful as negative constraints:**
> "the rates of catalytic oxidation reactions were independent of the presence of radical reaction
> inhibitors (KI, 10⁻² mol/l; HCOOK, 10⁻¹ mol/l) and activators (ammonia, ethylamine, or
> ethylenediamine, 10⁻¹ mol/l) for all of the three groups of compounds. In the presence of
> acetonitrile (up to 1 mol/l), the rates of oxidation of thiol compounds from the first group
> decreased several times (as well as the rates of oxidation of all the other thiols in weakly acidic
> and neutral media), whereas acetonitrile had no effect on the catalytic oxidation reactions of
> chelating aminothiols in an alkaline medium. The cyanide ion (2 × 10⁻³ g-ion/l) equally strongly
> inhibited the catalytic oxidation reactions of all the substrates."

**Cu²⁺ absence argument, verbatim (p. 379):** "the presence of Cu²⁺ (~10⁻⁴ mol/l) in aminothiol
solutions should be accompanied by intense absorption in the range 550–500 nm with **ε = 5000 l
mol⁻¹ cm⁻¹** [17, 18]. However, this is not the case." — ε value is **[C]** (Hanaki & Kamide 1974;
Hanaki 1980).

**Mechanism deferred, verbatim (p. 379):** "These mechanisms will be published elsewhere." — i.e.
this paper deliberately stops short of a rate-constant model.

## 10. DEFECTS AND CAVEATS TO CARRY FORWARD

1. **No rate constants.** Only initial rates w₀ with variable fractional orders. Eq. (3)'s k is never
   evaluated. The authors disclaim quantitative use in their own words (§2).
2. **No pKa values.** None anywhere.
3. **No temperature series → no Ea.** Single temperature, 20 °C.
4. **No error bars, no R², no per-point n.** Only a blanket "error of 5% or lower" on [O₂].
5. **No cysteine or glutathione rate data in any table.** Cysteine appears only in Figs. 1 and 3
   (figure-only) and in the group-2 prose ordering; GSH appears only in the prose ordering.
6. **The 250-fold HcySH/BuSH ratio has no tabulated support** and no stated conditions.
7. **TSN/ESN/DSN vs TSH/ESH/DSH transliteration split** between tables and text.
8. **Fig. 4's compound ((CH₃)HNCH₂CH₂SH) differs from Table 3's compound ((CH₃)₂NCH₂CH₂SH)** on the
   same page — do not merge them.
9. **Table 4 Block B has 8 rows and 9 order-entries**; the order-to-row mapping is not recoverable.
10. **Table 1 pH-13.0 prints two identical condition rows with different w₀** (3.35 and 1.68).
11. **Fig. 5 uses CuSO₄ (Cu²⁺)** where the rest of the paper uses [Cu(CH₃CN)₄]ClO₄ (Cu(I)).
12. Concentration units are **g-ion/l** for copper (a legacy Soviet unit numerically equal to mol/l
    for a monovalent ion) — carried through verbatim, not converted.

## NEW-PARAMETER TABLE (consolidated)

No rate constant, equilibrium constant, activation energy or pKa is available from this source. The
extractable parameters are reaction orders, initial rates, half-lives, one Brønsted coefficient, and
ordinal reactivity rankings.

| parameter | value | units (as printed) | conditions | anchor (table/page) | provenance |
|---|---|---|---|---|---|
| reaction order in [O₂], group 1 (TSH-type) | x = 1 throughout, except x = 1.00 held at all tabulated pH | dimensionless | 20 °C, phosphate–borate, pH 6–13 | Table 1, p. 374 / PDF p. 3; Eq. 4 p. 374 | [F] |
| reaction order in [RSH], group 1 | y decreases 1 → 0 with rising pH (tabulated: 0.73 @ pH 6; 0.50 @ 7; 0.25 @ 8; 0.03 @ 9; 0.11 @ 11; 0.63 @ 13) | dimensionless | as above | Table 1, p. 374 / PDF p. 3 | [F] |
| reaction order in [Cu⁺], group 1 | z rises fractional → 1 (tabulated: 0.22/0.55 @ 6; 0.25/0.75 @ 7; 0.59 @ 8; 0.83/1.00 @ 9; 0.95 @ 10; 0.95/1.00 @ 11; 1.00 @ 13) | dimensionless | as above | Table 1, p. 374 / PDF p. 3 | [F] |
| reaction order in [O₂], group 2 (ESH-type) | x: 1.00 @ pH 6–7.5; 0.90 @ 8; 0.75 @ 9; 0.35 @ 10; 0.10 @ 11; **0 @ 12 and 13** | dimensionless | 20 °C, phosphate–borate | Table 2, p. 376 / PDF p. 5; Eq. 5 p. 374 | [F] |
| reaction order in [RSH], group 2 | y: 0.60 @ 6; 0.81 @ 7; 0.48 @ 7.5; 0.42 @ 8; 0.36 @ 9; 0.35 @ 10; 0.34 @ 11; **0 @ 12 and 13** | dimensionless | as above | Table 2, p. 376 / PDF p. 5 | [F] |
| reaction order in [Cu⁺], group 2 | z: 0.25 @ 6; 0.30 @ 7; 0.52 @ 7.5; 0.60–1.00 @ 8; 1.00–1.25 @ 9; 1.50–1.75 @ 10; 1.90–2.00 @ 11; **2.00 @ 12**; 1.50 @ 13 | dimensionless | as above | Table 2, p. 376 / PDF p. 5 | [F] |
| reaction orders, group 3 (DSH-type) | x: 1.00 @ pH 6–9, 0.70 @ 10, 0 @ 11; y: 1.00 @ 6, 0.14 @ 7, 0.68 @ 8, 0.30 @ 9, 0 @ 10, 0.12 @ 11; z: 0.10 @ 6, 1.00 @ 7, 0.45 @ 8, 0.87 @ 9, 1.00 @ 10, 1.30 @ 11 | dimensionless | 20 °C, phosphate–borate | Table 3, p. 377 / PDF p. 6 | [F] |
| w₀ matrix, TSH | 24 values, 0.18 → 6.50 (× 10⁻⁷) | mol l⁻¹ s⁻¹ | 20 °C, [O₂]₀ 3 × 10⁻⁴ M, pH 6–13, [TSN]₀ 2.2–110 mM, [Cu⁺] 5–50 µM | Table 1, p. 374 / PDF p. 3 | [M] |
| w₀ matrix, ESH | 44 values, 0.15 → 596 (× 10⁻⁷) | mol l⁻¹ s⁻¹ | 20 °C, [O₂]₀ 3 × 10⁻⁴ M, pH 6–13, [ESN]₀ 1.0–100 mM, [Cu⁺] 1–500 µM | Table 2, p. 376 / PDF p. 5 | [M] |
| w₀ matrix, DSH | 18 values, 0.60 → 430 (× 10⁻⁷) | mol l⁻¹ s⁻¹ | 20 °C, [O₂]₀ 3 × 10⁻⁴ M, pH 6–11, [DSN]₀ 5–10 mM, [Cu⁺] 6.6–84.5 µM | Table 3, p. 377 / PDF p. 6 | [M] |
| w₀ matrix, TSH + ESH mixed | 18 values, 0.4 → 35.0 (× 10⁻⁷) | mol l⁻¹ s⁻¹ | 20 °C, pH 11, [O₂]₀ 3 × 10⁻⁴ M, [TSN]₀ 5–50 mM, [ESN]₀ 0–5 mM, [Cu⁺] 0.5–10 × 10⁻⁵ g-ion/l | Table 4, p. 378 / PDF p. 7 | [M] |
| τ½, O₂ binding to [Cu⁺(ESH)ₙ] | 2.90 | s | pH 8.5, [ESH] 5 × 10⁻² M, [Cu⁺] 1 × 10⁻⁴ g-ion/l | Fig. 6 caption, p. 379 / PDF p. 8 | [M] |
| τ½, O₂ binding to [Cu⁺(ESH)ₙ] | 0.78 | s | pH 11.5, same concentrations | Fig. 6 caption, p. 379 / PDF p. 8 | [M] |
| τ½, O₂ binding to [Cu⁺(PSH)ₙ] | 0.29 | s | pH 11.0, [PSH] 5 × 10⁻² M, [Cu⁺] 2 × 10⁻⁴ g-ion/l | Fig. 6 caption, p. 379 / PDF p. 8 | [M] |
| mixer dead time | ~30 | ms | jet-mixing polarographic setup | text p. 379 / PDF p. 8 | [M] |
| electrode response time | 2–3 | s | 2 µm gap, 10–15 µm Teflon membrane | Experimental p. 373 / PDF p. 2 | [M] |
| [O₂] measurement error | ≤ 5 | % | [O₂] = (2–100) × 10⁻⁵ mol/l | Experimental p. 373 / PDF p. 2 | [M] |
| Brønsted–Bjerrum coefficient | 1.02 | (dimensionless, with Z_aZ_b√I) | pH 10.5, µ = 0.1–0.5 mol/l, N-methylcysteamine | Eq. 6, p. 376 / PDF p. 5; Fig. 4, p. 377 / PDF p. 6 | [F] / [C] (Laidler, ref. 14) |
| HcySH : BuSH rate ratio | 250 | (dimensionless) | "at equal concentrations of Cu⁺"; **pH and concentrations NOT stated; no supporting table** | text p. 375 / PDF p. 4 | [M] (unsupported) |
| relative reactivity ordering, group 1 | HcySH > PSH > MfSH > PipSH > BSH > PentSH > TSH > EtSH > PrSH > BuSH | ordinal | 20 °C, Cu⁺-catalysed O₂ oxidation | text p. 374 / PDF p. 3 | [M] (ordinal only) |
| relative reactivity ordering, group 2 | NH₂CH(C₂H₅)CH₂SH ≈ NH₂CH(CH₃)CH₂SH > (CH₃)NHCH₂CH₂SH > ESH > **CySH** > TGA > **GSH** > NH₂C(C₂H₅)CH₂SH > NH₂CH₂CH₂CH₂NHCH₂CH₂SH | ordinal | as above | text p. 374 / PDF p. 3 | [M] (ordinal only) |
| relative reactivity ordering, group 3 | (C₂H₅)₂NCH₂CH₂SH > DSH > NH₂CH₂C(CH₃)₂SH > HOCH₂CH₂CH₂SH > HOCH₂CH₂SH | ordinal | as above | text p. 375 / PDF p. 4 | [M] (ordinal only) |
| ε, Cu²⁺–aminothiol, 500–550 nm | 5000 | l mol⁻¹ cm⁻¹ | used as a negative test for Cu²⁺ | text p. 379 / PDF p. 8 | [C] (refs. 17, 18) |
| thiol pKa | **NOT PRESENT IN THIS PAPER** | — | — | — | — |
| rate constant k | **NOT PRESENT IN THIS PAPER** | — | — | — | — |
| activation energy Ea | **NOT PRESENT IN THIS PAPER** | — | — | — | — |

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> These sources are not yet in `docs/reference/FIT_HOLDOUT_DECLARATION.md`. A declaration
> amendment is required before any wave may fit them. This section is a proposal only.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Everything in Tables 1–4 as *rate constants*** | **NOT FITTABLE — DECLARE INELIGIBLE** | — | The authors state in print that the rate law "can only serve as the qualitative characteristics" because the orders are variable and fractional. A w₀ with a fractional, pH-dependent order in three reactants is not convertible to a rate constant without inventing the missing order structure. Any parameter extracted from these tables would be a fabrication. |
| Reaction orders x, y, z vs pH for ESH (Table 2) — the transition x: 1 → 0 and z: 0.25 → 2.00 across pH 6 → 12 | **HOLD-OUT (structural/mechanistic, not numeric)** | pH | This is the paper's real content and it is a genuinely discriminating test: any Cu-catalysed thiol-oxidation submodel must reproduce the *switch* from an O₂-dependent, Cu-first-order regime below pH ~8 to an O₂-independent, Cu-second-order regime above pH ~11. Score as a qualitative regime-boundary prediction, never as a fitted number. |
| Reaction orders for TSH (Table 1) — x = 1 at every pH, z → 1, y → 0 | **HOLD-OUT (structural)** | thiol chelating ability | The orthogonal arm to the ESH row: a nonchelating thiol that never leaves O₂ first order. If the model reproduces ESH's switch but also switches for TSH, it is over-generalising. |
| The three ordinal reactivity series (§7) | **HOLD-OUT (ordinal only)** | thiol structure | 24 thiols ranked. Usable as a rank-correlation (Spearman) check on any structure–reactivity term the model grows. Explicitly **not** a magnitude check. Note CySH and GSH sit mid-series, which is itself informative for a food model. |
| τ½ = 2.90 / 0.78 / 0.29 s (Fig. 6 caption) | **HOLD-OUT** | pH (8.5 vs 11.5) and thiol identity (ESH vs PSH) | These are the only true time constants in the paper and are printed, not digitised. The ESH pair is a clean single-variable pH cut at fixed [ESH] and [Cu⁺]. |
| Ionic-strength result (µ = 0.1–0.5 has **no** effect on nonchelating thiols; accelerates chelating aminothiols in alkali; Brønsted coefficient 1.02) | **HOLD-OUT (sign/direction)** | ionic strength × chelating ability | A cheap, high-information directional test the trunk almost certainly does not encode. |
| Inhibitor/activator non-effects (KI, HCOOK, ammonia, ethylamine, ethylenediamine — all no effect; CN⁻ 2 × 10⁻³ g-ion/l — strong inhibition of all substrates; acetonitrile up to 1 M — several-fold suppression of group 1 but no effect on chelating aminothiols in alkali) | **FIT (as mechanism-selection constraints, not as parameters)** | — | These are qualitative gates on mechanism (rules out free-radical chain; requires Cu coordination). They set no numbers, so they cannot leak into a scored prediction. |
| The 250-fold HcySH : BuSH ratio | **DO NOT USE** | — | No supporting table, no stated pH, no stated concentrations. Unverifiable from this file. |
| Fig. 1 (cysteine time courses) and Fig. 3 (cysteine pH profile) | **DEFER — digitisation required** | pH | These are the only cysteine data in the paper and they are figure-only. If the wave needs a cysteine number, it must be digitised in a separate, explicitly [fig]-tagged pass with the render page recorded. I have not digitised them here. |

**Circularity risks flagged.**
(i) The reaction orders x, y, z in Tables 1–3 are **not independent of the w₀ values in the same
tables** — they were computed from those very rates via Eqs. (1)/(2). Fitting to w₀ and then scoring
against the orders (or vice versa) would be scoring the model against itself.
(ii) The group-1/2/3 classification is **defined by** the order behaviour, so "does the model put
thiol X in the right group?" and "does the model get thiol X's orders right?" are the same question
asked twice.
(iii) **No Ea and single temperature (20 °C).** Nothing here licenses any temperature transfer.
(iv) This is Cu(I)-catalysed autoxidation with a deliberately supplied [Cu(CH₃CN)₄]ClO₄ catalyst at
10⁻⁷–10⁻⁴ M and pH up to 13. Transferring any of it to a food matrix requires an explicit
free-copper-activity assumption the paper does not supply.
