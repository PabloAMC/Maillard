# Cremer 2000 — EXTRACTION (freeze-dried glucose 1.0 mol/kg + alanine, valine, isoleucine, leucine 50 mmol/kg each on microcrystalline cellulose, 0.1 mol/L citrate pH 5.0 before drying, aw 0.52, 80-110 C in sealed 22-mL headspace vials for 3-120 min; pseudo-zero-order formation of acetaldehyde, 2-methylpropanal, 2-methylbutanal, 3-methylbutanal by static headspace GC-FID with standard addition; four barriers; a glucose + leucine run at 90 C with Fru-Leu, leucine and A420; five plant powders at 90 C)
### The source of the "115 to 124 kJ/mol" Strecker barriers the corpus has carried second-hand since Jousse 2002: four amino acids in one low-moisture pot give barriers within 9 kJ/mol of each other, printed as text with a "5 %" standard deviation — and not one rate constant is tabulated.

**Source on disk:** `data/articles/rainercremer2000.pdf` (7 pp., 361,263 bytes, owner's download
2026-09-09; the file name carries the author's middle name — Dirk Rainer Cremer — as if it were a
surname; the citation is Cremer & Eichner 2000). Born-digital Elsevier PDF (Acrobat Distiller 3.01,
2000), clean text layer read from `scratchpad/articles/rainercremer2000.txt`; pages 39-42 rendered at
110 dpi to read the figure axes and printed bar labels. **The paper has no tables.** Every number below
is from the running text, the methods, or a label printed on a figure; the seven figures (plant powder
time courses at 90 C, Figs. 1-2; 60-min yields, Fig. 3; paprika free amino acids vs aldehydes, Fig. 4;
the glucose + leucine run, Fig. 5; the four-temperature isotherms, Fig. 6; the Arrhenius plot, Fig. 7)
are figure-only. No supplementary material (2000). Cremer's dissertation (Münster 1999) and the
companion paprika paper (Cremer & Eichner 2000, JAFC 48:2454) are cited and not on disk.

## 0. Identity

| field | value |
|---|---|
| Title | "The reaction kinetics for the formation of Strecker aldehydes in low moisture model systems and in plant powders" |
| Authors | Dirk Rainer Cremer*, Karl Eichner — Institut für Lebensmittelchemie der Universität Münster, Corrensstr. 45, 48149 Münster, Germany |
| Journal | Food Chemistry 71 (2000) 37-43; received 30 November 1999, accepted 29 February 2000 |
| DOI | 10.1016/S0308-8146(00)00122-9 (PII S0308-8146(00)00122-9 is printed; the DOI is the PII in Elsevier's scheme and is not itself printed on the page) |
| PDF file name | `data/articles/rainercremer2000.pdf` (misleading stem; not renamed) |
| On disk vs SI | full paper on disk; no SI exists |
| Naming | AA = acetaldehyde (from alanine); 2-MP = 2-methylpropanal (valine); 2-MB = 2-methylbutanal (isoleucine); 3-MB = 3-methylbutanal (leucine); G/L = glucose + leucine 20 : 1; G/AVIL = glucose + Ala + Val + Ile + Leu 20 : 1 : 1 : 1 : 1; Fru-Leu = the Amadori product fructose-leucine; ARP = Amadori rearrangement product |
| Who quotes it | Jousse 2002 (R9 Ea "115-124"; `jousse2002_extraction.md`), Balagiannis 2009 ("120 and 124 kJ/mol" for 3-MB and 2-MB — **swapped**, see Flag 4; `balagiannis2009_extraction.md`), Balagiannis 2015, Parker 2013 (ref 26), Huang 2017 |

## 1. Why it matters

The trunk's Strecker step `k_strecker` (rule R07) carries glycine only, and the amino-acid-identity wave
drafted in `results/validation/kinetic_core_b19_prereg_draft.md` needs per-amino-acid rates in water at
two or more temperatures; its row table lists this paper as "Cremer & Eichner 2000 (barriers 115 to
124)". That is exactly and only what the paper supplies: **four barriers, one per amino acid, from one
pot, at aw 0.52 on cellulose — not in water.** No rate constant is printed (Fig. 7 plots ln k on an
unlabelled scale). So the paper cannot give the wave a magnitude for any amino acid, and its barriers are
low-moisture barriers that the authors themselves place above the aqueous value (their comparison: 124
kJ/mol here vs Chan & Reineccius' 80.4 kJ/mol for 3-MB in water). What it does settle for the wave's
structure: (a) the four barriers agree within 9 kJ/mol (115, 115, 120, 124 with an estimated 5 % SD,
i.e. about +/- 6 each), so at aw 0.52 the barrier of the net Strecker cascade does not depend on the
amino acid among Ala, Val, Ile, Leu — support for the draft's "one barrier / one pH term shared across
amino acids" if it holds in water; (b) the glucose + leucine run shows the aldehyde appears only after the
Amadori product Fru-Leu has built to a steady state (about 40 min at 90 C) and is linear while Fru-Leu is
steady (30-120 min, r2 0.99), which is the mechanistic statement behind the trunk's structure: the
Strecker rate is set by the dicarbonyl supply from the Amadori pool, so a barrier measured on the
aldehyde in a sugar + amino-acid pot is the supply barrier. Compare `MARTINS_M4`: Amadori -> methylglyoxal
125.0 +/- 4.7 kJ/mol, Amadori -> 3-deoxyglucosone 97.0 +/- 1.7, Amadori -> 1-deoxyglucosone 107 +/- 7.3
(glucose + glycine, water, pH 6.8); and the pyrazine step's fed-dicarbonyl Strecker barriers 103.1 and
114.9 kJ/mol (`FROZEN_B18`). Cremer's 115-124 lands on Martins' Amadori-decomposition band, which is
what (b) predicts.

## 2. Methods as they matter to a model

- **Model G/AVIL (the Arrhenius pot).** "1.80 g glucose and 5.0 ml of each of the 0.1 M amino acid
  solutions (containing 44.5 mg Ala, 58.5 mg Val, 65.6 mg Leu, 65.6 mg Ile; sum of amino acids: 234
  mg), 20 ml buffer (pH 5.0; 0.56 g dry matter), 20 ml distilled water and (10.0 - 1.80 - 0.56 - 0.234 g
  =) 7.406 g microcrystalline cellulose were mixed, the pH corrected, using sodium hydroxide solution
  and/or hydrochloric acid, the mixture deep-frozen and freeze-dried." Stated result: **1.0 mol/kg
  glucose and 50 mmol/kg of each of the four amino acids** (per kg of dry matter). Check (mine): 1.80 g
  / 180.16 = 9.99 mmol glucose; 44.5 mg / 89.09 = 0.4995 mmol Ala; 58.5 / 117.15 = 0.4994 mmol Val; 65.6
  / 131.17 = 0.5001 mmol Leu and Ile; in 10.0 g dry matter -> 0.999 mol/kg and 50.0 mmol/kg. Buffer:
  **0.1 mol/L citrate, pH 5.0** (21.01 g citric acid monohydrate + 200 mL 1 M NaOH to 1000 mL; 96.4 : 3.6
  with 0.1 M NaOH), 20 mL = 2.0 mmol citrate -> **200 mmol/kg dry matter** (mine); the amino acids were
  dissolved "using as much hydrochloric acid as necessary", so an unstated amount of chloride is present;
  cellulose is 74 % of the dry mass. **The pH of the dry system is not a defined quantity; "pH 5.0" is
  the pH of the slurry before freeze-drying.**
- **Model G/L.** Same procedure, 1.0 mol/kg glucose + 50 mmol/kg leucine (20 : 1). Heated at 90 C only
  (Fig. 5).
- **Water activity.** 200 mg aliquots of the freeze-dried models (500 mg of plant powders) "adjusted to
  an aw-value of 0.52 by storing them for four days in a headspace vial placed in a desiccator over a
  saturated magnesium nitrate hexahydrate solution". Moisture content at aw 0.52 is not printed.
- **Heating = the headspace oven.** "The prepared sample aliquots were heated (thermostatted) for 3 to
  120 min in the headspace oven at 70 to 110 C" in sealed 22-mL Perkin-Elmer headspace vials (200 mg
  solid; the vial is essentially all headspace). Then the whole headspace is sampled: pressurisation
  0.8 min, injection 0.06 min, one injection per vial, so **each time point is a separate vial**. Fig. 6
  prints four isotherm labels, **80, 90, 100 and 110 C**, and Fig. 7's abscissa spans the same window;
  70 C appears only in the methods sentence (Flag 3). Time windows per isotherm are figure-only (Fig. 6
  shows about 90-150 min at 80 C down to about 10-15 min at 110 C).
- **Quantification.** Static headspace GC-FID, Perkin-Elmer HS-101 / 8410, 60 m x 0.32 mm x 1.0 µm
  Stabilwax, 40 C (5 min) to 70 C at 2 C/min (5 min); needle 120 C, transfer line 130 C. **Standard
  addition:** stock = 100 mg of each aldehyde in 100 mL diethylene glycol dimethyl ether (1.0 mg/mL
  each); 10, 20, 30 µL added to sample aliquots in the vials = 10, 20, 30 µg of each aldehyde per vial
  (mine) = for the 200 mg model 50 / 100 / 150 µg/g, i.e. 3-MB and 2-MB 0.58 / 1.16 / 1.74 mmol/kg, 2-MP
  0.69 / 1.39 / 2.08 mmol/kg, acetaldehyde 1.13 / 2.27 / 3.40 mmol/kg (mine; M = 86.13, 72.11, 44.05
  g/mol) — the spikes bracket the measured range (Fig. 6 ordinate runs to 1.8 mmol/kg). No internal
  standard; no response factors needed (standard addition in the matrix); no recovery, LOD or replicate
  count printed. Results in **mmol/kg** (of the aw-adjusted solid, presumably dry basis as charged).
- **Fru-Leu and leucine** (G/L run only): after the headspace analysis the vial contents extracted with
  1.0 mL water, centrifuged, Biotronic LC 5001 amino-acid analyser (Schräder & Eichner 1996 method),
  external calibration; Fru-Leu standard synthesised, 99 % pure by ion exchange. **Browning:** A420 of
  the same extract, 1.00 cm cell, diluted below 0.8 OD; Fig. 5 plots a "normalized absorption".
- **Plant powders.** Commercially dried cauliflower, spice paprika, asparagus, tomato, onion; 500 mg; aw
  0.52; 90 C; 15-120 min (Figs. 1-2), 60 min (Figs. 3-4). Free amino acids of paprika and tomato from
  Souci-Fachmann-Kraut 1989 and Cremer 1999, not measured here.
- **Kinetic analysis.** Zero-order slopes on the linear section of each isotherm; Arrhenius on four
  temperatures; "The r2 values of the linear regression of the isotherms in Fig. 6 ranged from 0.991 to
  0.999"; "The standard deviation of the activation energies was estimated to be at a 5 % level" (how is
  not said).

## 3. Tables re-typed

**The paper prints no tables.** What follows are the numbers printed in the text and as labels on the
figures, in the order they appear; the figures' plotted values are not typed.

### Printed in the text

| where | statement | number |
|---|---|---|
| Abstract, 3.3 | activation energies, model G/AVIL, aw 0.52, 80-110 C | **AA 115, 2-MP 115, 2-MB 120, 3-MB 124 kJ/mol** |
| 3.3 | SD of the Ea "estimated to be at a 5 % level" | about +/- 6 kJ/mol each (mine) |
| 3.3 | r2 of the zero-order isotherms of Fig. 6 | 0.991-0.999 |
| 3.2 | G/L at 90 C: zero-order section of 3-MB formation while Fru-Leu is at steady state | 30-120 min, r2 = 0.99 |
| 3.2 | G/L at 90 C: sigmoid 3-MB curve, inflection | about 180 min |
| 3.2 | G/L at 90 C: Fru-Leu reaches steady state | after about 40 min |
| 3.2 | G/L at 90 C: sum Leu + Fru-Leu + 3-MB | about 100 mol % of initial leucine until 120 min, "declined steeply" after |
| 3.3 | comparator: Chan & Reineccius 1994 (aqueous), Ea of 3-MB | 80.4 kJ/mol |
| 3.1 | plant powders: 2-MP, 2-MB, 3-MB zero order; acetaldehyde "parabolic" | — |
| 3.1 | acetaldehyde level "exceeds the level of the other aldehydes abundantly" though its rate is "comparable" | — |

### Digits printed as labels on figures (not axis readings; still figure-bound)

| figure | label | value | unit |
|---|---|---|---|
| Fig. 3 (60 min, 90 C, aw 0.52) | acetaldehyde produced, bar labels (bars exceed the axis) | cauliflower 7.5; spice paprika 2.3; asparagus 6.1; tomato 6.3; onion 1.6 | mmol/kg |
| Fig. 4 (spice paprika) | free alanine (literature value), bar label | 6.8 | mmol/kg |
| Fig. 4 | acetaldehyde produced after 60 min at 90 C, bar label | 2.3 | mmol/kg |
| Fig. 6 | isotherm temperature labels | 80, 90, 100, 110 | C |
| Fig. 7 | axis labels | "Ln k" (no unit) vs "1000/T, 1000/K" | — |

### Figures (all FIGURE-ONLY)

- Fig. 1: sum of 2-MP + 2-MB + 3-MB vs time (15-120 min) at 90 C in the five powders, ordinate
  "sum of aldehydes produced, mmol/kg", 0-1.6; straight lines through the origin region.
- Fig. 2: acetaldehyde vs time, same powders, 0-12 mmol/kg; concave-down curves for cauliflower,
  asparagus, tomato; near-linear low curves for paprika and onion.
- Fig. 3: the four aldehydes per powder at 60 min, 0-0.5 mmol/kg for the branched aldehydes.
- Fig. 4: paprika free Ala, Val, Ile, Leu vs the four aldehydes at 60 min.
- Fig. 5: G/L at 90 C, 0-360 min: Fru-Leu, Leu, their sum with 3-MB (left axis, mol % of initial
  leucine, 0-100); 3-MB and normalised A420 (right axis, mol %, 0-20).
- Fig. 6: G/AVIL isotherms at 80, 90, 100, 110 C, "aldehyde produced, mmol/kg", 0-1.8, four aldehydes
  per isotherm; the four lines within one isotherm lie close together (3-MB highest, acetaldehyde lowest
  at 100 and 110 C by eye).
- Fig. 7: Arrhenius lines for the four aldehydes; ordinate "Ln k" with no unit; the four lines are
  nearly parallel and nearly coincident. The scale of the ordinate is consistent with k expressed per
  gram and per second rather than per kilogram and per minute (my inference from the axis range against
  Fig. 6's slopes; recorded so that a later digitisation does not mis-assign the unit; no value typed).

## 4. Kinetic numbers the repository can use

Registry mapping: acetaldehyde -> `acetaldehyde`; 2-methylpropanal -> `2_methylpropanal`;
2-methylbutanal -> `2_methylbutanal`; 3-methylbutanal -> `3_methylbutanal`; **alanine, valine,
isoleucine, leucine, glucose, fructose-leucine (Amadori), citrate, cellulose -> not in
`data/keys/compounds.yml`.**

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| Strecker of leucine (net, supply-limited) | Ea, 3-methylbutanal formation | 124 (+/- ~6, "5 %") | kJ/mol | glucose 1.0 mol/kg + Ala, Val, Ile, Leu 50 mmol/kg each, citrate 200 mmol/kg (pH 5.0 before drying), cellulose carrier, aw 0.52, 80-110 C, sealed 22-mL vials, one vial per point | pseudo-zero order (slope of the linear section) | text 3.3, Fig. 7 | measured_barrier (low moisture; SD "estimated") |
| Strecker of isoleucine | Ea, 2-methylbutanal | 120 (+/- ~6) | kJ/mol | same | pseudo-zero order | text 3.3 | measured_barrier (low moisture) |
| Strecker of valine | Ea, 2-methylpropanal | 115 (+/- ~6) | kJ/mol | same | pseudo-zero order | text 3.3 | measured_barrier (low moisture) |
| Strecker of alanine | Ea, acetaldehyde | 115 (+/- ~6) | kJ/mol | same | pseudo-zero order | text 3.3 | measured_barrier (low moisture; acetaldehyde also has non-Maillard sources per the authors, in powders) |
| amino-acid identity, barrier | Ea(Leu) - Ea(Ala) | 9 | kJ/mol | same pot | — | text 3.3 | within_study_ratio (of barriers; inside the stated +/- 6 each) |
| amino-acid identity, rate | k(3-MB) : k(2-MB) : k(2-MP) : k(AA) at each T | — | — | same pot; lines near-coincident in Fig. 7 | pseudo-zero order | Fig. 7 | figure_only (the four rate constants are never printed) |
| zero-order rate constants, any aldehyde, any T | — | **not printed** | — | — | — | Fig. 7 only, unlabelled unit | figure_only |
| Amadori steady state (G/L, 90 C, aw 0.52) | time to steady Fru-Leu; linear 3-MB window; inflection of 3-MB | about 40 min; 30-120 min (r2 0.99); about 180 min | min | glucose 1.0 mol/kg + Leu 50 mmol/kg | — | text 3.2 | level_only (timescales printed as text) |
| leucine mass balance (G/L, 90 C) | Leu + Fru-Leu + 3-MB | about 100 mol % until 120 min | mol % of initial Leu | same | — | text 3.2 | level_only |
| Fru-Leu, Leu, 3-MB, A420 vs time (G/L, 90 C) | — | mol % | 0-360 min | — | Fig. 5 | figure_only |
| plant powders, 90 C, 60 min | acetaldehyde produced | 7.5 / 2.3 / 6.1 / 6.3 / 1.6 (cauliflower / paprika / asparagus / tomato / onion) | mmol/kg | aw 0.52, 500 mg | — | Fig. 3 bar labels | figure_only (printed digits, but on a figure; part of the acetaldehyde pre-exists in the powder per the authors) |
| plant powders, 90 C | branched aldehydes vs time; 60-min yields | — | mmol/kg | — | Figs. 1, 3, 4 | figure_only |
| comparator (not this paper's data) | Chan & Reineccius 1994, 3-MB in water | 80.4 | kJ/mol | aqueous, pH 6-8, 75-115 C | pseudo-zero order | text 3.3 | second-hand; see `chan1994b_extraction.md` (19.2 kcal/mol) |

**Can a second-order constant in water be derived?** No, on two independent grounds. (1) No rate
constant is printed: the ln k values exist only as points on Fig. 7 with an unlabelled ordinate. (2)
The system is a freeze-dried solid at aw 0.52; a zero-order rate in mmol per kg of dry solid per minute
has no volume and cannot be divided by a molarity. If a later wave digitises Fig. 7 (figure_only, +/- a
few tenths in ln k), the only licensed use is the **within-pot rate ratios among the four aldehydes at
one temperature** — the amino-acid identity factors the draft wants — under the assumption that the four
amino acids at equal 50 mmol/kg loading share one dicarbonyl pool, so that k(3-MB)/k(AA) =
k2(Leu)/k2(Ala) with the pool concentration cancelling. That ratio would transport to water only under
the further assumption that the relative reactivity of amino acids toward a dicarbonyl is the same in a
dry glass and in solution, which nothing on disk tests. The barriers transport as low-moisture
apparent barriers; the authors' own statement is that they exceed the aqueous value.

**Comparison with the trunk's constants.** Cremer's 115-124 kJ/mol vs the pyrazine step's
fed-dicarbonyl Strecker barriers 103.1 (glyoxal) and 114.9 (methylglyoxal) kJ/mol at pH 8, 100-120 C, and
vs Martins' Amadori-decomposition barriers 97-125 kJ/mol in water at pH 6.8. The agreement with the
Amadori band is what the authors' own reading of Fig. 5 implies (aldehyde rate = Fru-Leu decomposition
rate during the steady state); the agreement with the fed-dicarbonyl step is coincidental at this
precision (+/- 6 vs +/- 10). Against Chan 1994b's 80 kJ/mol (water, same aldehyde) the gap is 44
kJ/mol, and against Balagiannis 2009's 137 +/- 15 (liver extract, glucose -> intermediate) the gap is
-13; the ordering water < dry is the authors' claim (Hendel 1955), the ordering dry < liver extract is
not explained by moisture.

## 5. Flags

1. **No table of rate constants.** The four Arrhenius lines are the only record of k, and the ordinate
   carries no unit. Request from the authors (or Cremer's 1999 Münster dissertation, which will hold
   them): the zero-order k per aldehyde per temperature in mmol kg-1 min-1 with the fitted time windows,
   the number of vials per isotherm, and the SD computation behind "5 %".
2. **"5 %" SD is an estimate, not a fit statistic.** Four temperatures, one slope each; with r2 of
   0.991-0.999 on the isotherms the barrier's real uncertainty is the vial-to-vial scatter, unprinted.
   Carry +/- 6 kJ/mol as a floor.
3. **70 C in the methods, 80 C in the figures.** "70 to 110 C" (2.5) vs isotherm labels 80, 90, 100,
   110 C (Fig. 6) and an abscissa that stops near 80 C (Fig. 7). Either the 70 C series was run and not
   used, or the methods sentence is loose. Four temperatures at most.
4. **Balagiannis 2009 misquotes this paper**: it attributes "120 and 124 kJ/mol" to "3-methylbutanal
   and 2-methylbutanal, respectively"; the paper prints 2-MB 120, 3-MB 124. Jousse 2002 and Balagiannis
   2015 quote the range correctly.
5. **Low moisture, undefined pH.** aw 0.52 on 74 % cellulose; citrate 200 mmol/kg dry set to pH 5.0 in
   the slurry. The trunk has no water-activity term and no pH term below 6.8 on the sugar path; a rate
   from this pot cannot be placed on the trunk. The barriers are reported here as evidence about
   amino-acid (in)dependence, not as trunk values.
6. **Each time point is a separate vial and the whole headspace is injected.** The linear section is
   chosen after an induction period whose length is temperature-dependent (Fig. 6 starts its 110 C
   isotherm near 10 min, its 80 C isotherm near 90 min); the zero-order k is the steady-state slope, not
   an initial rate — consistent with the G/L reading (the aldehyde rate tracks the Amadori steady state)
   but it means the barrier includes the temperature dependence of the induction period only insofar as
   the window was chosen after it.
7. **Acetaldehyde in powders is partly pre-formed and partly from lipid / carotenoid oxidation** (the
   authors' own caveats, and its "parabolic" time course); its model-system barrier (115) is on a
   Maillard-only pot and is fine, but no acetaldehyde level from the powders is a Strecker yield.
8. **Standard addition detail.** Spikes of 10-30 µg per vial against 200 mg of solid: the three spike
   levels correspond to 0.6-1.7 mmol/kg for the C5 aldehydes, at or above the highest measured value
   (Fig. 6 tops at ~1.8 mmol/kg); the low end of the isotherms (0.1-0.2 mmol/kg) is an extrapolation of
   the addition line. Replicate count per point not printed.
9. **The G/L mass balance is a strong result printed only as prose.** "Sum ... remained approximately at
   the 100 mol % level" until 120 min, then fell: for the first two hours at 90 C every leucine that left
   is in Fru-Leu or 3-MB, so the Strecker yield per Amadori decomposed is high in a dry system with no
   competing water-mediated sinks — the opposite of Balagiannis 2009's liver extract (F_leu 2.3 %). The
   numbers behind this (Fig. 5) are figure-only.
10. **Registry gaps** (`data/keys/compounds.yml`): the four aldehydes are keyed; none of the four amino
    acids, glucose, the Amadori product or citrate is.
11. **Companion paper to fetch**: Cremer & Eichner 2000, JAFC 48:2454-2460 (paprika powder; Parker 2013
    ref 25) — same group, same method, may print the k values this paper omits.
