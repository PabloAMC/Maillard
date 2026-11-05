# Wietstock, Kunz & Methner 2016 (J. Agric. Food Chem. 64:8035) — EXTRACTION (file stem `wietstock2016`)
### *"Relevance of Oxygen for the Formation of Strecker Aldehydes during Beer Production and Storage"*
### Bottled **lager / Pilsner beer** (pH 4.41-4.44, 4.9 % ethanol, native Fe 48-53 µg/L) spiked with Leu / Ile / Phe (and **D5-Phe** 20 mg/L), glucose / xylose and Fe2+, bottled at **controlled O2** (total packaged oxygen 0.05-3.2 mg/L) and stored **20-28 °C for 12 weeks to 18 months**; plus a buffered model (pH 4.3, 5 % ethanol, Leu / Ile / Phe 25 mg/L, Fe2+ 100 µg/L, O2 0.49 mg/L, 6 months, 28 °C), mashing at 65 °C and wort heating to 100 °C under air vs N2. Strecker aldehydes by SAFE-GC/MS. **Strecker aldehydes rise linearly with total packaged O2** (sum of five: **4.5 ± 0.3 nmol per µmol O2**, i.e. ~220 mol O2 per mol aldehyde, in beer); **added Fe2+ has no effect, 10 mg/L EDTA suppresses D5-phenylacetaldehyde completely**. A free-amino-acid, metal / ROS route at storage temperature; **no Amadori product, no osone, no sugar effect**.

**Read 2026-11-05 (first reader) from `data/articles/wietstock2016.pdf`** (10 pp., ACS publisher PDF, clean text layer).
Prose and tables from `pdftotext -raw` (tables checked against 110 dpi renders); **Fig. 1 (p. 4) and Fig. 2 (p. 5)
digitised from 300 dpi crops** (scratchpad `sweep_img/w16f1-04.png`, `w16f2-05.png`; read against gridlines): read error
**±0.5 µg/L** (Fig. 1), **±0.05 µg/L** (Fig. 2). Fig. 4 and SI Fig. S1 from 110-130 dpi renders (±0.5 µg/L; ±0.05 mg/L).
Arithmetic: scratchpad `tools/sweep_fits.py` (maillard env, `-I`). Supporting information (`Wietstock2016_supplementary.pdf`,
11 pp.: Tables S1-S5, Figs. S1-S3) read into **sec. 8**.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | **Philip C. Wietstock** (corresponding), Thomas Kunz, Frank-Jürgen Methner |
| Affiliation | TU Berlin, Dept. of Food Technology and Food Chemistry, **Chair of Brewing Science**, Seestrasse 13, 13353 Berlin |
| Venue | *J. Agric. Food Chem.* **2016**, 64(42), 8035-8044 |
| DOI | **10.1021/acs.jafc.6b03502** (printed on every page) — matches the brief |
| Dates | received 4 Aug 2016; revised 4 Oct 2016; accepted / published 8 Oct 2016 |

**Laboratory.** Methner's brewing chair, TU Berlin = **`tu_berlin_methner`** in `data/keys/laboratories.yml` (strings
"TU Berlin Brewing (Jacob / Methner)"; Jacob 2018 / 2019). **Not** the Kanzler food-chemistry laboratory of the Bork
papers (different department and address). File name: the SI is `Wietstock2016_supplementary.pdf` (capital W);
aliased in `data/keys/pdf_dossier_aliases.yml` like the other capitalised SI stems.

## 1. What this paper answers for the model

| model need | what Wietstock 2016 gives | tag |
|---|---|---|
| Strecker aldehyde vs O2, quantitative | **Table S4**: five aldehydes + D5-PAL at 13 design points, bottling-atmosphere O2 < 0.01-20.99 % (TPO 0-3.2 mg/L), 12 weeks, 28 °C; linear in O2 (Table 5, S5) | [M] |
| yield per mol O2 | beer: **4.5 ± 0.3 nmol aldehyde (sum of 3-MB, 2-MB, methional, PAL, benzaldehyde) per µmol O2**; D5-PAL 0.37 ± 0.03 nmol/µmol; buffered model ≤ ~4.5 % (22 mol O2 per mol aldehyde, initial O2 only) | [Z] |
| Fe2+ effect | **none** (RSM: Fe2+ 0-500 µg/L insignificant, p > 0.05 for every aldehyde; Fig. 1 Leu/O2 vs Leu/O2/Fe2+ alike); **EDTA 10 mg/L blocks D5-PAL** and lowers the others | [M] |
| precursors | Leu, Ile, Phe (+ D5-Phe); glucose / xylose 131-468 mg/L (no effect on any Strecker aldehyde) | [M] |
| fed ARP, osone, dicarbonyl | **none** | — |
| heating with O2 at 100 °C | wort heated 20 → 100 °C (60 min + 5 min) after air vs N2 purge: aldehydes ×1.12-1.30 with air (Table 8) | [M] |

## 2. Conditions [M]

| experiment | matrix and charge | pH | T / time | vessel, headspace, O2 | metals / chelators |
|---|---|---|---|---|---|
| **A. RSM sugar / AA / Fe** (Table 1, S2) | lager beer 11.5 °P, 4.9 % v/v ethanol, 32 BU; + glucose / xylose 131.8-468.2 mg/L each, Leu / Ile / Phe 6.4-73.6 mg/L each, Fe2+ 0-100 µg/L (CCRD, 20 runs) | beer pH **4.41** (base beer; label, not re-measured) | **20 °C, 18 months**, dark | 500 mL bottles adjusted to 490 mL; overfoamed to expel headspace O2; crown cork | native Fe **48 µg/L** |
| **B. Leu / Fe / air** (Table 2, Fig. 1) | Pilsner 11.7 °P, 4.9 %, 28 BU; Leu 25 mg/L; or Leu / Ile / Phe 20 mg/L each; Fe2+ 500 µg/L | **4.44** | **28 °C, 30 weeks**, dark | overfoamed: **TPO 53.1 ± 3.3 µg/L**; air left in headspace: **2990.8 ± 177.3 µg/L** (n = 6 bottles) | native Fe **53 µg/L**; SO2 0.8-1.0 mg/L at spiking |
| **C. RSM O2 / Fe + D5-Phe** (Table 3, S4, Fig. 2) | Pilsner as B, pre-aged 20 weeks at 28 °C until SO2 < 0.1 mg/L; + **D5-Phe 20 mg/L**; Fe2+ 0-500 µg/L | 4.44 (base) | **28 °C, 12 weeks** | filled in a **N2 glovebox** whose O2 was set to < 0.01 / 3.07 / 10.5 / 17.92 / 20.99 %; TPO linear in that O2 (SI Fig. S1) | native Fe 53 µg/L |
| **D. EDTA / O2-barrier liner** (Table 6) | as C | 4.44 | 28 °C, 12 weeks | OB liner bottled at < 0.01 % O2; EDTA arm capped at 10.5 % O2 | **EDTA 10 mg/L** |
| **E. buffered model** (Fig. 4) | "0.1 mM acetate buffer" [!] (likely 0.1 M), **pH 4.3**, 5 % v/v ethanol, Leu / Ile / Phe 25 mg/L each (191 / 191 / 151 µM [Z]), Fe2+ 100 µg/L | 4.3 (set) | **28 °C, 6 months**, dark | N2-purged ≥ 30 min, 0.33 L glass bottles; **O2 490 ± 53 µg/L** (15 µM [Z]) after filling | Fe2+ 100 µg/L; **no anoxic or Fe-free control** [!] |
| **F. mashing** (Table 7) | 80 g Pilsner malt grist + 320 mL water, stirred | 5.50-5.51 (measured, cooled) | 65 °C, 90 min, then 78 °C | glovebox under N2 (< 0.01 % O2, deaerated water) vs air | — |
| **G. wort heating** (Table 8) | sweet wort 12.7 % extract, 250 mL in 250 mL Schott bottles | 5.57 (measured before and after, cooled) | 2 °C purge 2 h (air or N2, 100 mL/min), then 20 → 100 °C water bath in 60 min, +5 min at 100 °C | closed bottles; **dissolved O2 not measured** [!] | — |

## 3. Assays [M]

- **Aldehydes:** SAFE distillation, HRGC-MS (SIM ions in Table S1: m/z 58 for 2- / 3-methylbutanal, 104 methional, 91
  phenylacetaldehyde, 106 benzaldehyde), **pentanal as internal standard** (not isotope dilution); details in the
  authors' ref. 30 [C]. **Reporting limit "< 1.0 µg/L"** for every carbonyl (Tables 6-8, S4); no formal LOD / LOQ.
  - [!] Fig. 2 plots D5-PAL and D5-benzaldehyde values of 0.1-1.0 µg/L, below the tables' "< 1.0 µg/L"; the
    sub-1 µg/L points are below the stated reporting limit.
- **Replicate noise:** triplicate SDs in Tables 6-8: ±0.1-1.0 µg/L for the Strecker aldehydes (3-15 % relative);
  RSM pure error from 6 (A) or 5 (C) centre points: e.g. Table S4 centre (10.5 % O2, 250 µg/L Fe): 3-MB 11.1-13.7, PAL
  14.1-18.3, D5-PAL 2.8-3.7 µg/L (CV ~8-10 %).
- **TPO** by DIGOX 6.1 (MEBAK 2.28.1.1.2), n = 6 bottles; **Fe** by ICP-OES (matrix-matched); amino acids by EZ:Faast
  GC-FID; SO2 by CFA; radicals by ESR (spin trap POBN, EAP / T400).
- Recoveries: amino-acid and Fe additions 98.7 ± 8.3 % and 96.6 ± 0.8 %; D5-Phe and Fe in C > 95 %.

## 4. Data

### 4.1 Table S4 — RSM C, 12 weeks at 28 °C, D5-Phe 20 mg/L [M]; TPO from SI Fig. S1 [Z]

| O2 in atmosphere (%) | Fe2+ added (µg/L) | TPO (mg/L) [Z] | 3-MB | 2-MB | MET | PAL | BZALD | PAL-D5 | FUR |
|---|---|---|---|---|---|---|---|---|---|
| < 0.01 | 250 | ~0 | 7.2 | 4.4 | 2.7 | 10.0 | 1.8 | 1.6 | 474.0 |
| 3.07 | 427 | 0.46 | 9.3 | 5.6 | 4.0 | 11.4 | 2.9 | 1.9 | 487.1 |
| 3.07 | 73.2 | 0.46 | 9.4 | 5.5 | 3.3 | 10.9 | 2.6 | 1.9 | 444.4 |
| 10.5 | 0 | 1.59 | 12.1 | 7.0 | 4.5 | 15.6 | 5.4 | 3.3 | 476.3 |
| 10.5 (×5 centre) | 250 | 1.59 | 11.1-13.7 | 6.0-7.5 | 3.9-6.4 | 14.1-18.3 | 4.9-6.4 | 2.8-3.7 | 451-561 |
| 10.5 | 500 | 1.59 | 13.7 | 7.3 | 4.0 | 14.6 | 6.1 | 3.0 | 494.2 |
| 17.92 | 427 | 2.71 | 16.9 | 8.5 | 7.1 | 20.5 | 9.9 | 5.0 | 544.0 |
| 17.92 | 73.2 | 2.71 | 16.8 | 8.7 | 6.3 | 21.1 | 11.0 | 5.1 | 568.8 |
| 20.99 | 250 | 3.17 | 17.6 | 9.2 | 8.2 | 25.9 | 12.9 | 6.3 | 469.3 |
| fresh beer | — | — | 4.9 | 2.4 | < 1.0 | 1.1 | 1.8 | < 1.0 | 169.8 |

All in µg/L. D5-benzaldehyde < 1.0 everywhere. [!] The table's "fresh beer" row (MET < 1.0, PAL 1.1) differs from Table 6's
"reference (fresh)" (MET 2.8 ± 1.0, PAL 8.3 ± 1.5) for what the methods describe as the same base beer.

### 4.2 SI Fig. S1 — TPO vs O2 in the bottling atmosphere [F ±0.05 mg/L]

3.07 % → 0.47; 5.2 % → 0.77; 10.5 % → 1.58; 17.9 % → 2.71; 21 % → 3.18 (SD ≈ 0.15) mg/L; line through the origin, slope
**0.151 mg/L per % O2** (r² > 0.99 printed). 3.17 mg/L = **99 µM** O2 per litre of beer [Z].

### 4.3 Fig. 2 — D5-PAL and D5-benzaldehyde vs storage time, 28 °C [F ±0.05 µg/L]

| week | 0 | 1 | 2 | 3 | 5 | 7 | 12 |
|---|---|---|---|---|---|---|---|
| D5-PAL, 10.5 % O2 + 250 µg/L Fe2+ | 0.20 | 0.62 | 1.23 | 2.33 | 3.03 | 3.54 | 3.61 |
| D5-PAL, < 0.01 % O2, no Fe2+ | 0.09 | 0.31 | 0.51 | 0.71 | 0.82 | 1.05 | 1.30 |
| D5-benzaldehyde, 10.5 % O2 | 0.00 | 0.20 | 0.21 | 0.33 | 0.29 | 0.37 | 0.42 |

With O2 the curve is sigmoidal and plateaus by week 7 (O2 consumed); without, near linear [Q]. The 12-week values match
Table S4's 10.5 % rows (2.8-3.7) and < 0.01 % row (1.6) within ~0.3 µg/L.

### 4.4 Fig. 1 — experiment B, 30 weeks at 28 °C [F ±0.5 µg/L, week 30]

| | reference (overfoamed) | + Leu 25 mg/L | + Leu, air | + Leu, air, Fe2+ 500 µg/L |
|---|---|---|---|---|
| 3-MB | 8.7 | 12.0 | 25.0 | 23.0 |
| 2-MB | 6.4 | 6.2 | 9.9 | 9.5 |
| methional | 9.7 | 9.5 | 14.0 | 14.0 |
| PAL | 4.0 | 4.0 | 28.0 | 21.0 |
| benzaldehyde | 20.0 | 23.0 | 26.0 | 27.0 |

Text (vs reference): 3-MB +187 %, 2-MB +55 %, methional +44 % with air + Leu (agree with the figure: +187 / +55 / +44 %);
Leu alone +38 % 3-MB (agrees). [!] **Text PAL +64 % and benzaldehyde +583 % do not match Fig. 1** (figure: PAL +600 %,
benzaldehyde +30 %); the two values look swapped. Trust the figure.

### 4.5 Table 6 — 12 weeks, 28 °C [M, mean ± SD, n = 3]

| | fresh | aged reference | O2-barrier liner | EDTA 10 mg/L |
|---|---|---|---|---|
| 3-MB | 4.9 ± 0.7 | 7.5 ± 0.1 | 6.4 ± 0.2 | 5.5 ± 0.5 |
| 2-MB | 2.4 ± 0.4 | 4.7 ± 0.4 | 4.0 ± 0.2 | 2.9 ± 0.1 |
| methional | 2.8 ± 1.0 | 4.9 ± 1.0 | 3.1 ± 0.3 | 4.2 ± 0.2 |
| PAL | 8.3 ± 1.5 | 14.2 ± 0.1 | 12.4 ± 0.9 | 10.0 ± 0.7 |
| benzaldehyde | 1.2 ± 0.3 | 2.0 ± 0.1 | 1.8 ± 0.2 | 1.4 ± 0.2 |
| **PAL-D5** | < 1.0 | **1.6 ± 0.3** | < 1.0 | **< 1.0** |

[!] The aged reference's O2 at bottling is not stated (the EDTA arm used 10.5 %); its PAL-D5 1.6 sits at Table S4's
< 0.01 % level, not its 10.5 % level (2.8-3.7).

### 4.6 Fig. 4 — buffered model, 6 months, 28 °C [F ±0.5 µg/L; SD]

3-MB 7.5 (2.5), 2-MB 7.7 (2.8), PAL 22.6 (0.7), benzaldehyde 33.8 (1.5) µg/L (text: "7.5 ... to 33.8 µg/L" [M]).
Methional not measured (no Met added).

### 4.7 Tables 7 / 8 — wort [M, mean ± SD, n = 3]

Mashing, air vs N2 (65 °C, 90 min): 3-MB 75.1 vs 26.8; 2-MB 29.7 vs 20.3; methional 39.7 vs 13.2; PAL 47.0 vs 14.4;
benzaldehyde 7.1 vs 3.9 µg/L (×1.5-3.3); free amino N 210 vs 219 mg/L; pH 5.51 vs 5.50.
Heating to boil after air vs N2 purge: 3-MB 51.4 vs 39.5; 2-MB 21.4 vs 16.7; methional 29.1 vs 26.0; PAL 39.9 vs 34.7;
benzaldehyde 7.8 vs 6.4 µg/L (×1.12-1.30); pH 5.57 in both, unchanged by heating.

### 4.8 Arithmetic [Z] (`sweep_fits.py`)

1. **Yield per O2 in beer (Table S4 + Fig. S1, 13 points, 12 weeks, 28 °C):** linear slopes in nmol aldehyde per µmol
   TPO: 3-MB 1.23 ± 0.09, 2-MB 0.53 ± 0.05, methional 0.48 ± 0.08, PAL 1.24 ± 0.13, benzaldehyde 1.05 ± 0.08, **sum 4.5 ±
   0.3** (→ **~220 mol O2 consumed per mol Strecker aldehyde made**); D5-PAL 0.37 ± 0.03 nmol/µmol (from 117 µM D5-Phe:
   ≤ 0.04 % converted). The rest of the O2 goes to SO2, polyphenols, hop acids, ethanol (hydroxyethyl radical) [Q / P].
2. **Buffered model:** sum 0.68 µM aldehydes vs 15.3 µM initial O2 → **4.5 % per O2, 22 mol O2 per mol aldehyde**.
   Direction of error: crown-cork O2 ingress over 6 months raises the true O2 (more O2 per aldehyde); aldehyde losses
   (oxidation to acids) lower it; no anoxic arm, so part of the aldehyde may not need O2. Not bounded either way.
3. **Wort heating:** excess aldehydes with air ≈ 0.28 µM (sum of five); dissolved O2 not measured. If the air-purged
   wort held its air saturation at 2 °C (~0.4 mM O2 in water, handbook value; less in 12.7 % extract), that is
   ≥ ~1000 mol O2 per mol extra aldehyde — an order of magnitude only.
4. **Fe2+:** native Fe 48-53 µg/L (~0.9 µM) already saturates the O2-limited channel; the authors' reading [Q]. EDTA
   10 mg/L (~34 µM) suppresses it.

## 5. What it can and cannot test in the model

**Can:** show that **free amino acids yield Strecker aldehydes from O2 without a dicarbonyl or Amadori step** in an
acidic, metal-containing matrix (D5-Phe → D5-PAL linear in O2, EDTA-suppressible, sugar-independent) — a third
O2-dependent form beside the M3 case's routes A (dicarbonyl) and B (Amadori oxidative). Give an order of magnitude for
its O2 cost (tens to hundreds of mol O2 per mol aldehyde).

**Cannot:**
- **Transfer to model pots.** Beer is a reductant-rich matrix (SO2, polyphenols, iso-α-acids, 4.9 % ethanol) at pH 4.4
  and 20-28 °C over weeks to months; the O2 stoichiometry is set by what else consumes O2. The engine's M3 pots are
  phosphate-buffered sugar / ARP / amino-acid solutions at ~100 °C, pH ~7, for minutes to hours. Neither the rate
  (weeks at 28 °C) nor the yield per O2 carries over; the direction (O2 + metals → aldehyde from free amino acid)
  might.
- **Test route B or R2 directly**: no ARP, no osone, no dicarbonyl measured; the authors cite direct ARP oxidation as a
  possible beer contributor without testing it [C].
- The buffered model has no anoxic, no Fe-free and no EDTA arm, and its buffer strength is misprinted or very weak.
- Wort heating (the only ~100 °C step) did not measure O2.

## 6. Conflicts and agreements with the repo's records

1. **Agrees in direction with Hofmann & Schieberle 2000b Table 2** (fed ARP-Phe: PA air / argon 9.2, air + Cu higher;
   M3 case D3 / F17#6) and Voigt 2010 (O2 raises a dicarbonyl-route aldehyde): O2 and transition metals raise Strecker
   aldehydes. Wietstock adds that **a metal-catalysed free-amino-acid route exists** and is EDTA-suppressible; so an
   air / argon or + Cu ratio is not by itself a route-B marker (the M3 case already says this for Voigt).
2. **[!] No sugar effect** on any Strecker aldehyde in beer at 20 °C / 18 months (Table 4): not a conflict with the
   engine (no record at 20 °C), but a caution that a dicarbonyl route is not visible at storage temperature in beer.
3. **[!] Internal:** text vs Fig. 1 percentages for PAL / benzaldehyde (sec. 4.4); fresh-beer rows of Table 6 vs Table
   S4 (sec. 4.1); aged-reference O2 not stated (sec. 4.5); sub-reporting-limit points in Fig. 2 (sec. 3).
4. No ledger or fit record cites this paper; no measured value is contradicted.

## 7. Proposed ledger entry (strecker / M3 context; NOT written — first reader only, second reader pending)

```yaml
  - id: wietstock2016_tableS4_strecker_aldehydes_vs_packaged_oxygen_beer_28c
    dossier: wietstock2016_extraction.md
    paper: "Wietstock, Kunz & Methner 2016 J Agric Food Chem 64:8035"
    lab: "TU Berlin Brewing (Jacob / Methner)"
    system: "Pilsner beer pre-aged to SO2 < 0.1 mg/L, + D5-phenylalanine 20 mg/L, Fe2+ 0-500 ug/L, bottled at 0-21 % O2 (TPO 0-3.2 mg/L)"
    medium: "beer pH 4.44, 4.9 % ethanol, native Fe 53 ug/L; 490 mL in 500 mL bottles, crown cork"
    temperature_c: [28]
    ph: {initial: [4.44], measured_during: false, final: null}
    time_min: [120960]
    species: [3-methylbutanal, 2-methylbutanal, methional, phenylacetaldehyde, benzaldehyde, phenylacetaldehyde-d5]
    kind: endpoint
    location: "SI Table S4 p. 6-7; TPO from SI Fig. S1 (dossier sec. 4.1-4.2)"
    status: digitised_in_dossier
    quant_class: calibrated   # pentanal internal standard, not an isotopologue
    simulable: false
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "sample bottles were finally prepared in a glovebox under nitrogen atmosphere"
    flags: [beer_matrix, storage_temperature, not_simulable, reporting_limit_1ugL, fe_insignificant, edta_suppresses]
    note: "Sum of five Strecker aldehydes 4.5 +- 0.3 nmol per umol packaged O2 (~220 mol O2 per mol aldehyde); D5-PAL 0.37. Order-of-magnitude context for M3's O2 forms only."
```

**Provenance tiers (first reader only):** Table S4, Tables 6-8 `direct_measurement` (printed, internal standard,
not SIDA), confidence medium; Figs. 1, 2, 4, S1 values `direct_measurement` digitised; every yield per O2 in sec. 4.8 is
[Z] `derived`, beer-matrix-specific, not a band for model pots.

## 8. Supplementary material (read 2026-11-05, first reader)

**Source.** `data/articles/Wietstock2016_supplementary.pdf` (11 pp.). Identified by its first page: this paper's title,
the three authors and the brewing-chair address. [!] The PDF's metadata title is "Antioxidant properties of hops during
wort boiling" (a leftover from another manuscript); the content is this paper's SI. Text by `pdftotext -layout`;
Fig. S1 from a 130 dpi render.

| item | content | tag |
|---|---|---|
| Table S1 | GC-MS retention times and quantifier ions; **pentanal is the internal standard** | [M] |
| Table S2 | RSM A raw data (20 runs, 18 months at 20 °C): 3-MB 6.2-11.0, 2-MB 3.2-8.5, MET 4.3-7.0, PAL 6.5-12.0, BZALD 1.8-3.4, FUR 399-742 µg/L; centre-point spread 3-MB 6.7-8.0, PAL 7.2-11.0 | [M] |
| Table S3 | RSM A model statistics (R² 0.20-0.87) | [M] |
| **Table S4** | RSM C raw data (sec. 4.1) | [M] |
| Table S5 | RSM C model statistics (R² 0.82-0.93; lack of fit n.s.) | [M] |
| **Fig. S1** | TPO vs O2 in the bottling atmosphere (sec. 4.2) | [F] |
| Fig. S2 | ESR signal vs time of the fresh beers of experiment B: air raises the radical signal at once; Fe2+ raises T400 | [Q] |
| Fig. S3 | time courses of ten carbonyls, 10.5 % O2 + 250 µg/L Fe2+ vs < 0.01 % O2, no Fe2+, 12 weeks, 28 °C: Strecker aldehydes sigmoidal with O2, near linear without; FUR, NON, NIC, heptanal, octanal alike in both | [Q] |

## 9. Open questions it bears on

**M3's O2 budget.** Rival R2 (route A only, with an O2-made osone) needs **≥ ~250-700 µmol O2 per mmol ARP** in air (case
D3: ≥ ~0.5-0.7 mol free osone per mol ARP within 120 min).

| item looked for | result |
|---|---|
| quantitative Strecker aldehyde formation vs O2 (headspace or dissolved) | **FOUND** — beer, 28 °C, 12 weeks: linear in total packaged O2 0-3.2 mg/L (0-99 µM); sum of five aldehydes **4.5 ± 0.3 nmol per µmol O2**, D5-PAL 0.37 ± 0.03 (Table S4 + Fig. S1). Model buffer: ~4.5 % per initial O2. Wort heating to 100 °C: ×1.12-1.30 with air, O2 not measured |
| Fe2+ effects | **FOUND (null)** — added Fe2+ 0-500 µg/L insignificant in both RSM designs and in Fig. 1; native Fe ~50 µg/L; **EDTA 10 mg/L suppresses D5-PAL completely** and lowers 3-MB / 2-MB / PAL by ~25-40 % vs the aged reference (Table 6) |
| amino acids and precursors | **FOUND** — Leu, Ile, Phe (6-74 mg/L each), D5-Phe 20 mg/L; glucose / xylose 132-468 mg/L (no effect on Strecker aldehydes); no ARP, no dicarbonyl, no Met added (methional from beer's own Met) |
| a rate or yield per mol O2 usable as a band | **NOT FOUND for M3's pots.** The yields exist (4.5 nmol / µmol O2 in beer; ≤ ~45 nmol / µmol in the model) but belong to a different route (metal / ROS oxidation of free amino acid), a reductant-rich matrix, pH 4.3-4.4 and 20-28 °C over weeks. They are not the O2 cost of making an osone from an ARP at 100 °C, which is what R2's 250-700 µmol O2 / mmol ARP needs |
| ARP or osone with O2 | **NOT FOUND** |

**Transfer to model pots:** direction only. **Effect on M3's O2 forms [Z]:** nothing changes R2's requirement or its
numbers. The paper does document a **third O2-dependent form**, free amino acid + O2 (+ trace metal) → Strecker
aldehyde, EDTA-suppressible and sugar-independent, which is in neither route A nor route B of the case. Whether it
matters in a 100 °C phosphate pot (Hofmann 2000b's air + Cu arm is the nearest test) is an owner call; no rate here
transfers.
