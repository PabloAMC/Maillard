# Wang & Arntfield 2015 — EXTRACTION (mixed C6-C8 aldehydes / ketones on salt-extracted pea and canola isolates, 1 % w/v, pH 8, and 95 C heating 0-60 min; ITEX headspace GC/MS + DSC)
### The heat paper of the Manitoba series: the only time course of aldehyde binding to a plant isolate on disk, printed as numbers for canola but not for pea.

**Source on disk:** `data/articles/wang2015.pdf` (owner's download, 2026-09-08). Read from the
scratchpad text layer (`wang2015.txt`, clean); pypdf confirms one table (Table 1, DSC of pea
isolate, p. 5) and four figures. Repo status before this dossier: `data/lit/binding_constants.yml`
cites this paper (Food Hydrocoll. 43:410-417) as thesis Chapter 4 under source
`wang_2015_umanitoba_thesis` but carries no record from it; its pea records come from the FRI 2015
salts/pH chapter (Table 5.1) and the unpublished hydrolysis chapter (Table 7.2). The engine's
matrix layer (`results/validation/matrix_sites_prereg.md` §4) concludes that covalent aldehyde
binding "over cooking times is a percent-level effect"; this paper's canola numbers are the direct
empirical check on that statement and disagree with it (see §1, §5).

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Binding of selected volatile flavour mixture to salt-extracted canola and pea proteins and effect of heat treatment on flavour binding" |
| Authors | Kun Wang, Susan D. Arntfield (University of Manitoba) |
| Venue | Food Hydrocolloids 43 (2015) 410-417; received 17 Jan 2014, accepted 19 Jun 2014 |
| DOI | 10.1016/j.foodhyd.2014.06.011 (PII S0268-005X(14)00237-9) |
| Naming | CPIs / PPIs = salt-extracted canola / pea protein isolates (same preparations as Wang 2014); "1-hexanal" = hexanal; "homological" = homologous; "heterological" = mixed-class (hexanal + 2-hexanone). "Control" in Figs 1-2 = the single-flavour, unheated system. |
| Companion | Wang & Arntfield 2014, Food Chem. 157:364-372 (`wang2014_extraction.md`): identical isolates, buffer, loading, headspace method. |
| Compound registry | hexanal -> `hexanal`; heptanal -> `heptanal`; octanal, 2-hexanone, 2-heptanone, 2-octanone -> not in registry. |

## 1. Why it matters

Need (a). The engine binds hexanal to protein amine sites with a pseudo-first-order rate bracket
whose anchors are dairy-protein adduct studies at ambient to mild temperature, and reports
percent-level binding over cooking times. This paper heats a 1 % plant isolate with hexanal at 95 C
and follows headspace depletion at 0, 0.5, 1, 2, 5, 10, 20, 40, 60 min. For the canola isolate it
prints the numbers: hexanal bound goes from 14.34 % to 44.25 % in the first 30 s and to 66.85 % by
10 min, then stays; the ketone 2-octanone rises 17.99 -> 28.4 % by 2 min and falls back to 17.48 %
by 60 min (release on aggregation), which the authors read as aldehyde binding being irreversible
(covalent) and ketone binding reversible. For the pea isolate only the 95 C / 30 min end point is
measured, in mixtures, and the absolute values are FIGURE-ONLY; the text prints the increments
(+14.35 to +18.25 percentage points for the three aldehydes). This is the closest thing on disk to a
rate of aldehyde uptake by a plant globulin under heat; it is canola, it is depletion not adduct
yield, and it is a within-study time series, so under the repository's rule it may inform a rate
bracket only as a within-study shape, never as a pea constant. Table 1 (DSC) is fully printed and
gives the pea isolate's denaturation temperature (89.98 C) and enthalpy with and without flavours.
Need (b) is not touched.

## 2. Methods as they matter to a model

- **Isolates.** PPIs from commercial yellow pea flour by the 0.3 M NaCl / dilution / dialysis route
  of Wang 2014; CPIs by PMM from canola meal (Burcon AL018 this time). Protein (Dumas, N x 5.7):
  **CPIs 87.32 %, PPIs 82.68 %** — the same batches or the same numbers as Wang 2014.
- **Buffer, loading, flavour.** 2 % w/v isolate in 0.01 M potassium phosphate pH 8, ultrasonicated
  20 min. Single-flavour vials: 1 mL 2 % protein + 0.5 mL buffer + 0.5 mL of 1000 ppm stock
  (0.1 mL/100 mL, volumetric) -> **1 % w/v isolate (10 g powder/L = 8.27 g protein/L for PPIs),
  250 ppm v/v flavour (~204 mg/L; hexanal 2.03 mM, heptanal 1.79 mM, octanal 1.60 mM; ketones the
  same to two figures)**, 2 mL liquid in a 20-mL crimp vial. Mixtures: 1500 ppm stocks
  (0.15 mL/100 mL), 1/3 mL of each of three + 1 mL 2 % protein -> 250 ppm of each, **750 ppm total**
  (three-flavour) or 500 ppm total (hexanal + 2-hexanone, by the same construction). Flavour added
  last.
- **Equilibration.** 30 C, 125 rpm, 3 h ("Preliminary testing found that 3 h was adequate to reach
  equilibrium"). Duplicate vials, each sampled once.
- **Heat treatment (verbatim):** "For the heat treatment, reaction vials were placed into an
  Isotemp Water Bath 2320 (Fisher Scientific, Marietta, OH, USA) at 95 C and heated for the
  specified heating time. The control was not heated." Two experiments: (1) CPIs + hexanal and
  CPIs + 2-octanone, single flavour, heating time **0, 0.5, 1, 2, 5, 10, 20, 40, 60 min**; (2) one
  treatment, **95 C / 30 min**, for the three mixture systems on both CPIs and PPIs. Flavour was
  present during heating ("heating CPIs-flavour mixtures at 95 C"). Whether the 3-h 30 C
  equilibration preceded or followed the heating step is not stated (flag 3). Cited denaturation
  temperatures: CPIs ~89 C, PPIs ~86 C, "lower than the 95 C that was used".
- **Headspace sampling and GC/MS.** As Wang 2014: "After mixing, samples were incubated and shaken
  for 14 min at 40 C and 1 mL of sample headspace was aspirated into the GC injector port by a
  CombiPal autosampler unit with PAL Itex-2 (In-Tube-Extraction) absorber attachment (CTC Analytics
  AG, Switzerland) after one absorption cycle." Varian CP-3800 / 320-MS, VF-5ms 30 m x 0.2 mm,
  He 4 mL/min, 25 C/min to 265 C, hold 3 min; EI 70 eV, m/z 25-250. No internal standard, no
  calibration.
- **Quantification sentence (verbatim):** "Binding percentage of flavours was determined from the
  difference between the peak areas of flavoured samples in the absent and presence of proteins such
  that: Binding % = (1 - Peak area with protein added / Peak area without protein added) x 100%".
  For heated vials the protein-free reference is presumably heated the same way, but this is not
  stated (flag 4).
- **DSC (PPIs only).** 10 % w/v PPIs in 0.3 M NaCl with 250 ppm of each flavour (single, or each
  of the mixture), 1 h rotary shaking; 10-15 uL hermetic Tzero pans; 30 -> 120 C at 10 C/min;
  duplicates. Enthalpy unit printed as "J/K" in Table 1 and text; Wang 2014 printed J/g for the
  same measurement (flag 5).
- **Statistics.** Tukey p < 0.05; superscripts in Table 1 and the figures.

## 3. Tables re-typed

### Table 1. "Effect of flavour binding on the thermal properties of salt-extracted pea proteins."

| Flavour types | Carbon number | Denaturation temperature (Td, C) | Enthalpy of denaturation (delta-H, J/K) |
|---|---|---|---|
| Control (a) | — | 89.98 +/- 2.18 a | 14.53 +/- 0.14 a |
| Aldehyde flavours (b) | C6 | 88.56 +/- 0.10 a | 11.27 +/- 0.08 cd |
| | C7 | 91.12 +/- 4.23 a | 10.28 +/- 0.40 de |
| | C8 | 86.59 +/- 0.02 a | 8.87 +/- 0.16 e |
| Ketone flavours (b) | C6 | 92.91 +/- 0.99 a | 12.85 +/- 0.35 b |
| | C7 | 91.61 +/- 2.84 a | 11.93 +/- 0.05 bc |
| | C8 | 89.76 +/- 0.19 a | 11.00 +/- 0.62 cd |
| Mixture of three aldehydes (c) | — | 86.22 +/- 0.09 a | 7.64 +/- 0.47 f |
| Mixture of three ketones | — | 88.18 +/- 0.47 a | 11.35 +/- 0.04 c |
| Mixture of hexanal and 2-hexanone | — | 88.38 +/- 0.70 a | 10.76 +/- 0.08 cd |

Footnotes as printed: a-f, column values with the same letter not significantly different
(P < 0.05). (a) Control was not added with flavour compounds. (b) In the simple protein-flavour
system, only one flavour was added at 250 ppm; C6 indicates the respective aldehyde or ketone with
six carbons. (c) In the binary or ternary systems, each flavour was added at 250 ppm. Conditions:
10 % w/v PPIs in 0.3 M NaCl (not the binding buffer). Re-summed check: the text's quoted values
(7.64; 11.27, 10.28, 8.87; 11.35; 12.85, 11.00; 10.76) all match the table.

### Numbers printed only in the running text (Figs 1-4 are otherwise FIGURE-ONLY)

| system | quantity as printed | where |
|---|---|---|
| CPIs + hexanal, 95 C time course | "In the first 30 s, the percentage of hexanal bound to CPIs significantly increased by about 32% from 14.34 to 44.25% ... With further heating, the retention of hexanal continuously increased up to 66.85% over the first 10 min and then remained constant." | §3.2.1 (Fig. 3) |
| CPIs + 2-octanone, 95 C time course | "In the first 2 min of heating, the retention of 2-octanone increased from 17.99 to 28.4% ... binding dropped from 28.40 to 17.48% in the following 58 min." | §3.2.2 (Fig. 4) |
| CPIs, ternary aldehydes, unheated | hexanal binding "significantly increased" vs single (no number); "a much lower decrease in the flavour retention for octanal (16.57%) compared to heptanal (30.45%) in the ternary system" | §3.1.2 (Fig. 1b) |
| PPIs, ternary aldehydes, unheated | "the binding of all aldehydes was significantly increased by 16.36-20.40% for PPIs compared with the control" | §3.1.2 (Fig. 1d) |
| CPIs, ternary aldehydes, 95 C / 30 min vs unheated ternary | "heating enhanced binding of hexanal by 18.35% and heptanal and octanal by about 11.06 and 8.33%, respectively" | §3.2.3 (Fig. 1b) |
| PPIs, ternary aldehydes, 95 C / 30 min vs unheated ternary | "aldehyde flavour retention was also higher for PPIs, with increases between 14.35 and 18.25%" | §3.2.3 (Fig. 1d) |
| CPIs, hexanal + 2-hexanone, 95 C / 30 min | "binding of hexanal to CPIs increased from 24.85% to 33.44% but this was accompanied by a dramatic decrease in the binding of 2-hexanone from 14.9% to 3.68% (Fig. 2A). A similar pattern was followed by PPIs (Fig. 2B)." | §3.2.5 |
| Ketones, ternary, heated | CPIs: competition promoted (2-octanone favoured); PPIs: "the binding of ketones to PPIs was nearly unaffected" | §3.2.4 (Fig. 1a, c) |
| Ternary ketones, unheated | 2-octanone > 2-heptanone > 2-hexanone on both proteins; "a significant increase in the binding of 2-octanone by PPIs" in the mixture; others comparable to control | §3.1.1 |
| Cited (Damodaran & Kinsella 1981, soy) | nonanal 1094 M^-1 > 2-nonanone 930 M^-1 > 5-nonanone 541 M^-1 | §3.1.2 |
| Cited (Kuhn 2008, WPI 80 C) | 2-octanone: continuous decrease 1-80 min, 26 % remained bound | §3.2.2 |
| Cited (Ng 1989; Gkionakis 2006) | heated faba PMM (95 C / 15 min) and soy (60 C): 35.6 % and 20 % higher binding of vanillin / lactones | §3.2.3 |

The "%" increments are read as percentage points (e.g. 24.85 -> 33.44 is +8.59 points); the paper
does not say so, and its own "about 32%" for 14.34 -> 44.25 (a 29.91-point rise) fits neither a
points nor a relative reading exactly (flag 2).

## 4. Numbers the repository can use

| quantity | value | unit | conditions | source location | evidence class | `binding_constants.yml` fit |
|---|---|---|---|---|---|---|
| Hexanal % bound to CPIs, unheated (t = 0 of the heat series) | 14.34 | % | 10 g CPIs powder/L (8.73 g protein/L), 250 ppm v/v (2.03 mM), 0.01 M K-phosphate pH 8, 30 C / 3 h, headspace 40 C | §3.2.1 text | measured (single point printed in text) | `percent_bound_at_conditions`, protein_source canola isolate (no canola matrix exists in the repo; the record type itself fits) |
| Hexanal % bound to CPIs after 95 C for 0.5 min | 44.25 | % | as above + 95 C water bath, flavour present | §3.2.1 | measured | same type, with `protein_state: heated 95 C 0.5 min`; **canola** |
| Hexanal % bound to CPIs after 95 C for 10-60 min (plateau) | 66.85 | % | as above | §3.2.1 | measured | same; canola |
| 2-Octanone % bound to CPIs: 0 / 2 / 60 min at 95 C | 17.99 / 28.4 / 17.48 | % | as above | §3.2.2 | measured | same; canola; ketone |
| Hexanal % bound to CPIs, hexanal + 2-hexanone mixture, unheated -> 95 C / 30 min | 24.85 -> 33.44 | % | 250 ppm each (500 ppm total), otherwise as above | §3.2.5 | measured | same type; competitive system, canola |
| 2-Hexanone % bound to CPIs, same mixture, unheated -> heated | 14.9 -> 3.68 | % | as above | §3.2.5 | measured | same; canola |
| **PPIs**: hexanal / heptanal / octanal % bound, single flavour, unheated | FIGURE-ONLY | % | 10 g PPIs powder/L (8.27 g protein/L), 250 ppm v/v, pH 8, 30 C / 3 h | Fig. 1d "control" bars | figure_only | would be the wanted pea `percent_bound_at_conditions` rows; **cannot be filled** |
| PPIs: same, ternary mixture, unheated and 95 C / 30 min | FIGURE-ONLY | % | 250 ppm each | Fig. 1d | figure_only | as above |
| PPIs: increment of aldehyde binding, ternary mixture vs single flavour (unheated) | +16.36 to +20.40 | percentage points (reading, flag 2) | as above | §3.1.2 | measured (difference only) | no record type: an increment without its base |
| PPIs: increment of aldehyde binding on heating 95 C / 30 min (ternary) | +14.35 to +18.25 | percentage points (reading) | as above | §3.2.3 | measured (difference only) | `denaturation_effect_evidence` is the closest existing bucket (Barallat-Perez 2023 sits there qualitatively); this would be its first numeric entry, as an increment, pea, mixture |
| CPIs: heating increments, ternary aldehydes | hexanal +18.35, heptanal +11.06, octanal +8.33 | percentage points (reading) | as above | §3.2.3 | measured (difference) | as above, canola |
| CPIs: competition decrements, ternary aldehydes unheated | octanal -16.57, heptanal -30.45 | percentage points (reading) | as above | §3.1.2 | measured (difference) | none |
| PPIs Td (control) | 89.98 +/- 2.18 | C | 10 % w/v PPIs, 0.3 M NaCl, 10 C/min | Table 1 | measured | none; a matrix datum (denaturation onset for a pea isolate) |
| PPIs delta-H (control) | 14.53 +/- 0.14 | "J/K" as printed (J/g by the 2014 paper) | as above | Table 1 | measured | none |
| PPIs delta-H with 250 ppm hexanal / heptanal / octanal | 11.27 / 10.28 / 8.87 | as above | as above | Table 1 | measured | none; flavour-induced partial unfolding, 22-39 % enthalpy loss |
| PPIs delta-H with three aldehydes together | 7.64 +/- 0.47 | as above | 250 ppm each | Table 1 | measured | none |
| CPIs / PPIs cited Td | ~89 / ~86 | C | Uruakpa & Arntfield 2005; Sun & Arntfield 2012 | §2.6 | level_only, secondary | none |

### Within-study shape the rate bracket may consult (canola, depletion, not a pea constant)

From the printed CPIs-hexanal series (14.34 % at 0, 44.25 % at 0.5 min, 66.85 % from 10 min):
the rise above baseline reached at 30 s is (44.25 - 14.34) / (66.85 - 14.34) = 29.91 / 52.51 =
0.570 of the eventual rise, so the half-time of the approach to plateau is **under 30 s** at 95 C
(if single-exponential, k = -ln(1 - 0.570) / 0.5 min = **1.7 min^-1**, t1/2 ~ 25 s). Read on the
free fraction instead (0.8566 -> 0.5575 in 0.5 min), the first-order constant is 0.86 min^-1. Either
way it is orders of magnitude faster than the engine's ambient adduct brackets extrapolated to 95 C
would give (prereg §4: 0.1 % bound for hexanal in 1 % BLG at 160 C / 30 min), and it plateaus at
two-thirds bound rather than going to completion, which points to a finite site pool (or a
partition equilibrium of the unfolded protein), not to a simple first-order sink. These are
derived numbers on canola and are for the flag, not for a record.

## 5. Flags

1. **All pea percent-bound values are FIGURE-ONLY** (Fig. 1c, 1d, 2B). The text prints pea only as
   increments; a pea `percent_bound_at_conditions` record cannot be written from this paper.
2. **"%" increments are ambiguous** between percentage points and relative change; the paper's own
   "about 32%" for a 29.91-point rise does not settle it. Read as points here, flagged.
3. **Order of heating and equilibration not stated**; also not stated whether the 3-h 30 C step was
   applied at all to the heated vials, or whether the "0 min" point of Figs 3-4 is the 3-h
   equilibrated sample. The 14.34 % baseline is consistent with an unheated equilibrated control.
4. **Protein-free reference under heat not described.** If the no-protein vial was not heated
   identically, part of the "increase in binding" on heating is loss of hexanal to the heated
   headspace/septum or to heat-driven side reactions. Aldol self-condensation of hexanal at 95 C /
   pH 8 in the absence of protein is not excluded by anything printed.
5. **Enthalpy unit printed "J/K"** in Table 1 and text; the 2014 paper prints J/g for the same
   DSC protocol. Treated as J/g of sample; the Td column is unaffected.
6. **Depletion ≠ covalent adduct.** The irreversibility argument (hexanal keeps rising while
   2-octanone is released as the protein aggregates) is the authors' inference; no adduct was
   measured. The repository's covalent brackets and this depletion series measure different things,
   and the gap between them (§4 shape) should be recorded as a tension, not resolved by fitting.
7. **Canola is not pea.** Wang 2014 found canola binds aldehydes more than pea and pea binds
   ketones more than canola; the canola heat series must not be transferred to pea. The pea heat
   information in this paper is one point (95 C / 30 min), in a 750-ppm mixture, as an increment.
8. **"ppm" is volumetric** (0.1 mL / 100 mL stock), so 250 ppm is ~204 mg/L, not 250 mg/L — the
   same caveat as Wang 2014 flag 3, relevant to the `flavor_concentration_mg_per_L: 250.0` fields of
   the thesis-derived records already in `binding_constants.yml`.
9. **DSC medium differs from binding medium** (10 % protein in 0.3 M NaCl vs 1 % in 0.01 M
   phosphate); the enthalpy data support the unfolding story qualitatively but at 10x the protein
   and 30x the salt.
10. **Gao 2020 cites "Wang & Arntfield 2015"** for "PPI extracted at alkaline condition showed a
    higher binding capacity to aldehydes"; that finding is Wang 2014 (PPIa > PPIs), and the 2015
    paper Gao lists is the FRI salts/pH paper, not this one. Noted so the three 2014-2015 Wang
    papers are not conflated in the repo.
