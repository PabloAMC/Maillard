# Novotny, Cejpek & Velisek 2008 — EXTRACTION
### *"Formation of Carboxylic Acids during Degradation of Monosaccharides"*
### Glucose, fructose, arabinose, DL-glyceraldehyde, dihydroxyacetone (0.1 M), NO amine, refluxed 1 h in (1) K2S2O8 / water, (2) K2S2O8 / 0.3 M NaOH, (3) 0.05 M NaOH; formic, acetic, propionic acid and 24 hydroxy acids + 12 lactones by GC/MS

**Read 2026-10-31 from `data/articles/novotny2008.pdf`** (15 pp., Czech J. Food Sci. 26:117-131; born-digital text
layer). No ledger entry and no fit rows by decision (`data/lit/fit_rows/ph_acid_davidek_acids.yml` header and
`data/lit/source_verification/ph_acid_davidek_acids.yml`, last block: every condition is persulfate or NaOH without
amine).

**Legend.** **[V]** verbatim · **[M]** printed number · **[Q]** qualitative statement · **[C]** cited · **[Z]** derived
by me · **[!]** flag.

## 1. What this paper answers for the model

Question (B7a-2, `results/validation/kinetic_core_dpo_partition_prereg.md` sec. 11): does sugar -> acetic acid need
an amine; pH dependence; labelling; yields.

| need | what Novotny 2008 gives | tag |
|---|---|---|
| (a) amine-free | **All runs are amine-free.** In 0.05 M NaOH (1 h reflux, no oxidant) acetic acid is the major acid: **Glc 11.48, Fru 12.06, Ara 12.45, Dha 7.81, Gla 10.64 % (n/n)**; formic acid 5.32 / 5.91 / 5.65 / 6.47 / 6.69 % | [M] |
| (b) pH / oxidant | acetic acid is low with persulfate: K2S2O8 / water (acidifying) 0.52-2.57 %, K2S2O8 / 0.3 M NaOH 1.21-2.56 %; formic acid is high with persulfate (7.6-18.0 % and 8.8-20.6 %). No pH measured anywhere [!] | [M] |
| (c) labelling | **none** (no isotopes); mechanisms argued from product patterns and from Davidek 2006a / b [C] | — |
| (d) conditions / method | sugar 5 mmol in 50 mL (0.1 M) of water + K2S2O8 5 mmol, 0.3 M NaOH + K2S2O8 5 mmol, or 0.05 M NaOH; reflux (~100 °C [Z]) 1 h; 3 independent determinations; RSD < 20 %. Volatile acids: pH 3 (HCl), 4 x 20 mL ether, heptadecane IS, 5-point external calibration (sodium formate, sodium acetate, propionic acid), GC/MS (CP-WAX 52 CB), ions m/z 46 / 43 / 74 | [M] |

**Answer.** An amine-free, oxidant-free alkaline pot turns ~11-12 mol % of a hexose or pentose (arabinose) into acetic
acid in 1 h at reflux, so **acetic acid formation from sugars does not need an amine** in alkali (confidence ~95 %
in the observation; the route is not identified: no labels, no deoxyosone measured). It cannot price the amine
effect at near-neutral pH: there is no amine arm and no buffered near-neutral condition, so it neither supports nor
refutes B7a-2's amine-assisted cleavage. [Z] The base is exhausted during the run: total acids in 0.05 M NaOH reach
86.7-89.7 % (n/n) for Glc / Fru / Ara (~0.087 M acid against 0.05 M NaOH), so these pots end below neutral and the
"alkaline" yield is an integral over a falling pH. The pentose acetic acid ceiling it offers (~12 mol %, arabinose)
is well below Davidek 2008's xylose + glycine yields (16 / 28 mol % at pH 6 / 8, 120 °C, 4 h), but the conditions
differ in temperature, time, buffer and pH trajectory, so the comparison is not an amine effect estimate.

## 2. Identity

| field | value |
|---|---|
| Authors | Ondrej Novotny, Karel Cejpek, Jan Velisek (Department of Food Chemistry and Analysis, Institute of Chemical Technology, Prague) |
| Venue | *Czech J. Food Sci.* **2008**, 26 (2), 117-131 (as printed) |
| DOI | **10.17221/2465-CJFS** (not printed in the PDF; Crossref 2026-10-31 match on title, journal, volume 26, year 2008). [!] Crossref gives pages 113-131; the PDF prints 117-131 on every page |
| Dates | received 30 Jan 2007; accepted after corrections 2 Apr 2007 |
| Content | Tables 1-3 (yields % n/n per model); Figs. 1-2 (alpha- / beta-cleavage schemes after Davidek 2006b), 3-7 (proposed routes to formaldehyde, acetaldehyde, glyoxylic, propionic acid) |

## 3. Yields [M] (% n/n of sugar, averages of 3; RSD < 20 %)

| acid | model | Glc | Fru | Ara | Dha | Gla |
|---|---|---|---|---|---|---|
| formic | K2S2O8 / H2O | 9.93 | 7.59 | 18.04 | 10.66 | 11.63 |
| acetic | K2S2O8 / H2O | 0.79 | 1.68 | 0.52 | 2.39 | 2.57 |
| formic | K2S2O8 / 0.3 M NaOH | 20.64 | 13.50 | 15.64 | 8.76 | 10.11 |
| acetic | K2S2O8 / 0.3 M NaOH | 1.53 | 1.81 | 1.21 | 1.97 | 2.56 |
| formic | 0.05 M NaOH | 5.32 | 5.91 | 5.65 | 6.47 | 6.69 |
| acetic | 0.05 M NaOH | 11.48 | 12.06 | 12.45 | 7.81 | 10.64 |
| propionic | 0.05 M NaOH | 0.16 | 0.19 | 0.20 | 0.05 | 0.14 |
| lactic | 0.05 M NaOH | 7.07 | 7.64 | 5.08 | 10.40 | 7.50 |
| glycollic | 0.05 M NaOH | 14.51 | 16.30 | 26.02 | 10.39 | 12.59 |
| total acids + lactones | 0.05 M NaOH | 89.7 | 88.7 | 86.7 | 54.9 | 77.9 |

Totals in the other models: K2S2O8 / H2O 34.2 / 12.9 / 43.7 / 35.0 / 22.9; K2S2O8 / 0.3 M NaOH 98.3 / 92.6 / 124.4 /
116.7 / 129.6 % (n/n). Formaldehyde and acetaldehyde identified (OPD benzimidazoles), not quantified.

[Q] Acetic acid "correlated well with the amounts of acetol and methylglyoxal" (their 2007 paper), and the authors
propose acetaldehyde oxidation / Cannizzaro and radical recombination as extra routes besides the Davidek
beta-cleavage [C].

## 4. Flags

- [!] The Fig. 2 caption labels A1 / A2 as "hydrolytic alpha-dicarbonyl cleavage" and B1 / B2 as "oxidative
  beta-dicarbonyl cleavage", the reverse of the text on p. 119 and of Davidek 2006b, from which it is adopted.
- [!] No pH, no temperature beyond "reflux", no time course, no blank without sugar.
- Hydroxy acids for which no standard existed were semi-quantified against "structurally similar" acids.
- Fit table cross-check (2026-10-31): the Novotny numbers quoted in
  `data/lit/source_verification/ph_acid_davidek_acids.yml` (Tables 1-3, formic / acetic, all five sugars) match the PDF.
