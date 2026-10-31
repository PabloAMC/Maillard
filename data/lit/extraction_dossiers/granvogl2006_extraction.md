# Granvogl, Bugan & Schieberle 2006 — EXTRACTION
### *"Formation of Amines and Aldehydes from Parent Amino Acids during Thermal Processing of Cocoa and Model Systems: New Insights into Pathways of the Strecker Reaction"*
### Leu / Ile / Phe / **Met** (not Val) + methylglyoxal (2-oxopropanal), equimolar 50 µmol, **170 °C**, 2.5-30 min: model A on silica gel (10 % water), model B in 0.1 M phosphate pH 7.0 (5 mmol/L each). Strecker aldehyde by SIDA, "Strecker amine" by dansyl-LC-MS/MS. Plus unfermented / fermented / roasted cocoa (Sulawesi, Ghana).

**Read 2026-10-31 from `data/articles/granvogl2006.pdf`** (10 pp., JAFC 54:1730-1739; born-digital, clean text layer;
Tables 6-8 checked against a 110 dpi render of p. 1736, every cell matches). Figures 1, 3, 4, 5 are reaction schemes,
no numbers. No row of this paper is keyed anywhere in the repository yet.

**Legend.** **[V]** verbatim (short quote) · **[M]** number printed in this paper · **[Q]** qualitative statement by the
authors · **[C]** cited from another paper · **[Z]** derived by me, never printed · **[!]** flag.

**[!] The brief and `docs/guides/MECHANISM_REVIEWS_2026-10-31.md` sec. 3 describe this paper as "Ile / Leu / Phe / Val +
MGO".** The model systems use **L-Ile, L-Leu, L-Met, L-Phe** (Methods, p. 1732; Tables 6-9). Valine appears only in the
cocoa part (2-methylpropylamine, Table 2; free Val, Table 5). **Methylpropanal is measured nowhere in the paper**, neither in
cocoa nor in a model.

## 1. What this paper answers for the model

Context: the forced-input Strecker layer (`results/validation/kinetic_core_strecker_forced_input_prereg.md`) has one
Strecker flux per donor with two exits (aldehyde, acid), no amine exit, f_Leu declared 1 +- 0.3 dex, and covers 98 - 100 °C,
pH 7 only; MECHANISM_REVIEWS sec. 1 finding 7 says a fed-donor yield is a branching share, not a rate.

| need | what Granvogl 2006 gives | tag |
|---|---|---|
| (a) aldehyde time courses from MGO + each amino acid | Tables 6-9: 5 time points (2.5, 5, 10, 20, 30 min) x 2 models, mmol per mol amino acid, RSD of triplicates. Aqueous model B at 30 min: **3-methylbutanal 84.6, 2-methylbutanal 122.1, phenylacetaldehyde 72.8, methional 22.5 mmol/mol** (peak 25.3 at 20 min) | [M] |
| (b) the corresponding amines | Aqueous model B at 30 min: **3-methylbutylamine 0.26, 2-methylbutylamine 0.13, 2-phenylethylamine (PEA) 0.78, 3-(methylthio)propylamine 0.72 mmol/mol**. Low-water model A: 3.44 / 2.43 / 6.13 / 4.71 | [M] |
| (c) aldehyde : amine ratio | Model B, 30 min: **Leu 325, Ile 940, Phe 93, Met 31**; model A: Leu 52, Ile 165, Phe 4.3 (printed "4:1"), Met 2.7. The ratio FALLS with time: x1.7 - 30 from 2.5 to 30 min in model B (at 2.5 min: Leu 900, Ile 1590, Phe 450, Met 910), x10 - 300 in model A. The amine appears with a lag, the aldehyde without one. The 2.5-min amines (0.01 - 0.04) sit near quantitation, so the early ratios are loose | [M] / [Z] |
| (d) amine as a share of aldehyde + amine (the would-be amine exit) | model B, 30 min: **Leu 0.3 %, Ile 0.1 %, Phe 1.1 %, Met 3.0 %**; model A: 1.9 / 0.6 / 19 / 27 % | [Z] from [M] |
| (e) amino-acid identity, MGO donor, water, pH 7, 170 °C | aldehyde at 30 min relative to Phe: **Ile 1.68, Leu 1.16, Phe 1, Met 0.31**; at 10 min 1.57 / 1.12 / 1 / 0.44. Ile / Leu 1.44 (30 min), 1.41 (10 min) | [Z] from [M] |
| (f) conditions | 50 µmol amino acid + 50 µmol MGO; model B in 10 mL 0.1 M phosphate pH 7.0 (5 mmol/L each), "closed glass vessels", 170 °C; model A homogenised with 3 g silica gel 60 "containing 10% of water", 170 °C, closed vessels. Atmosphere NOT controlled or stated; heating medium and heat-up time NOT stated; pH set before heating, not re-measured | [M] / [!] |
| (g) quantification | aldehydes: SIDA ([2H2]-2-methylbutanal, [2H2]-3-methylbutanal, [2H3]-methional, [13C2]-phenylacetaldehyde), GC-MS (CI, methanol). PEA: SIDA ([2H4]-PEA). Other amines: internal standard [2H9]-n-butylamine with response factors 0.54 - 4.90 (**not SIDA**). Triplicates, RSD 0.5 - 37 % | [M] |
| (h) Strecker acids, MGO course, O2 | **none measured.** No acid (phenylacetic / 3-methylbutanoic), no MGO time course, no air / argon contrast | [!] |
| (i) amine -> aldehyde back-route | PEA + MGO (model A conditions): **~0.8 mol % phenylacetaldehyde and 7.1 mol % benzaldehyde** ("data not shown") | [M] (text only) |

**Answer.**

1. **In a water pot with methylglyoxal at pH 7 the Strecker amine is a negligible exit.** It is <= 1.1 % of the aldehyde
   for Leu / Ile / Phe and 3 % for Met at 30 min, 170 °C. Confidence in the numbers: ~0.9 (SIDA aldehyde, LC-MS/MS amine,
   triplicates). Confidence that the amine exit is also < 2 % of the Strecker flux at 98 - 100 °C, pH 7 (the layer's scope):
   ~0.75. That would be a transfer across 70 K, and the amine's lag suggests a slower secondary route. **The prereg's
   two-exit structure (aldehyde, acid; no amine) is supported for alpha-dicarbonyl donors in water.** The amine matters in low
   water (model A: Phe 19 %, Met 27 %) and in roasted cocoa (PEA exceeds PA ~50-120x molar, sec. 4.4). The lipid-enal class
   gives the same contrast: Zamora finds the amine is the product without O2 (MECHANISM_REVIEWS strecker review). An amine
   exit belongs to a low-moisture / matrix arm, not to step 2.
2. **Identity factor at the donor level: Leu / Phe ~1.1 - 1.2 with MGO in water at 170 °C** [Z]. log10 1.16 = 0.06 dex, inside
   the prereg's declared f_Leu = 1 +- 0.3 dex. Hofmann 2000 Table 2 gives Leu / Phe 2.1 (aldehyde) / 1.94 (aldehyde + acid)
   in the GLUCOSE pot at 98 °C (`hofmann2000_extraction.md`). The two lie 0.25 dex apart: a temperature effect, a supply effect
   (glucose pot vs fed MGO), or the missing acids (Granvogl measures the aldehyde only). Not separable here. Ile / Leu 1.4
   (model B) agrees with Balagiannis 2009's 1.6. Met is low and decays after 20 min (methional loss), so f_Met cannot be read
   from these yields without a methional sink.
3. **No rate constant, no branching share of the MGO flux.** MGO is never measured, and neither are the acids. The
   aldehyde-only yields bound s_MGO from below only for the aldehyde exit. Model B Ile approaches a plateau of ~12 % by 20 - 30 min
   (117.6 -> 122.1) and Phe (63.0 -> 72.8) is still rising. At 98 °C, 100 mmol/L, Hofmann 2000 Table 3 reaches ~4 % PA + PAA
   at 500 min. [Z] A crude lifetime argument puts MGO at tau ~10 - 15 min at 170 °C and [AA] = 5 mmol/L. That gives
   kappa_MGO ~ 1/(tau [AA]) ~ 1.5e-2 L/(mmol min), against ~5e-4 at 98 °C from Hofmann's 15 - 25 min at 100 mmol/L, i.e.
   x30 over 98 -> 170 °C, an apparent Ea of ~65 kJ/mol. **Order of magnitude only** (P ~0.3 that a proper fit lands within
   +-30 kJ/mol of this): unknown heat-up, uncontrolled O2, amine-free MGO sinks that matter more at 5 mmol/L, and different
   amino acids. Usable as a hold-out for a later temperature law (roadmap step 3), not as a step-2 input.
4. **Decarboxylation pathways (Fig. 3, proposed [Q]).** Amino acid + alpha-dicarbonyl gives an imine. Decarboxylation gives
   two tautomers: hydrolysis of tautomer 1 yields the aldehyde (+ aminoketone); hydrolysis of tautomer 2 yields the amine
   (+ the dicarbonyl back). Two further reasons for the low amine in model pots are proposed [Q]:
   - excess MGO traps the amine as a Schiff base, which converts it to aldehyde (Fig. 4 route a; measured 0.8 mol %) or to
     benzaldehyde + formaldehyde (route b, aza-vinylogous retro-aldol; 7.1 mol %);
   - aldimines resist hydrolysis (Rizzi 1974 [C]).

   For cocoa they propose an Amadori-direct oxidative decarboxylation via 3-oxazolines (Fig. 5) that can release either
   aldehyde or amine. This is hypothesis only, not measured.

   The lag is in the data [Z from M]. From 2.5 to 30 min in model B, the amine grows x6.5 - 72 (Ile 6.5, Leu 13, Phe 39,
   Met 72) while the aldehyde grows x2.5 - 8 (Met 2.5, Ile 3.8, Leu 4.7, Phe 8.1).
   That fits an amine made on a slower path, or one that survives only once MGO is depleted (the authors' trapping argument).
   It does not fit a fixed tautomer split.
5. **pH, buffer, air / argon: not varied.** One pH (7.0, 0.1 M phosphate), no atmosphere control, so nothing here on the
   pH factor g_D or on omega.

## 2. Identity

| field | value as printed |
|---|---|
| Authors | Michael Granvogl (TU München), Susanne Bugan (DFA), Peter Schieberle* (TU München and Deutsche Forschungsanstalt für Lebensmittelchemie, Garching) |
| Title | as in the header |
| Venue | *J. Agric. Food Chem.* **2006**, 54 (5), 1730-1739 |
| DOI | **10.1021/jf0525939** (printed p. 1730 and as the article id "JF0525939" on p. 1739) |
| Dates | received 19 Oct 2005; revised 29 Dec 2005; accepted 29 Dec 2005; published on web 31 Jan 2006 |
| Content | Table 1 (LC-MS/MS transitions, response factors), 2 (amines in cocoa), 3 (amine odour thresholds), 4 (aldehydes in cocoa), 5 (free amino acids in cocoa), **6-9 (aldehyde and amine vs time, models A / B)**; Figs. 1 (enzymic decarboxylation), 2 (chromatograms), 3 (Strecker scheme, aldehyde vs amine), 4 (PEA + MGO -> PA / benzaldehyde), 5 (Amadori -> 3-oxazoline hypothesis) |
| Companion | ref 10, Schieberle, Koehler & Granvogl 2005, Adv. Exp. Med. Biol. 561:205 (3-aminopropionamide from Asn; possibly `data/articles/schieberle2005.pdf`, not checked) |

## 3. Methods [M]

| item | as printed |
|---|---|
| Model A (low water) | L-Ile, L-Leu, L-Met or L-Phe (50 µmol) "singly reacted with 2-oxopropanal (50 µmol) by homogenizing both reactants with silica gel (3 g, silica gel 60, 0.063-0.200 mm) ... (containing 10% of water) and heating at 170 °C in closed glass vessels". Then water 15 mL + labelled amine + labelled aldehyde, 10 min equilibration, ultrasonication, centrifugation. [Z] ~0.3 g water per 3 g gel, so ~0.17 mol/kg water per reactant if all dissolved; a_w not stated |
| Model B (aqueous) | same four amino acids (50 µmol each) "singly dissolved in phosphate buffer (0.1 mol/L; pH 7.0; 10 mL)" and heated at 170 °C in closed glass vessels after adding MGO (50 µmol), i.e. **5 mmol/L each**. Internal standards added after cooling |
| Times | 2.5, 5, 10, 20, 30 min (both models) |
| Not stated | heating device / medium, heat-up time, headspace volume, atmosphere (no purge mentioned for the model pots; the argon purge in Methods belongs to the dansylation of reference amines), final pH, MGO purity / assay |
| Aldehydes (models) | 5 mL aliquot on Extrelut NT 20, eluted pentane / ether 9:1, concentrated, GC-MS. SIDA with [2H2]-2-methylbutanal, [2H2]-3-methylbutanal, [2H3]-3-(methylthio)propanal, [13C2]-phenylacetaldehyde |
| Amines (models) | 5 mL aliquot + NaHCO3 + dansyl chloride (3 h, dark), CH2Cl2 extraction, LC-MS/MS (ESI+, SRM; Table 1). [2H4]-PEA for PEA (SIDA); [2H9]-n-butylamine with response factors for the others: 2-methylpropyl 0.54, 2-methylbutyl 0.80, 3-methylbutyl 0.90, PEA 0.92 (vs [2H4]-PEA), 3-(methylthio)propyl 4.90. LoQ 0.2 µg/kg (text p. 1733); "limit of detection ... (0.1 µg/kg)" (p. 1734) |
| Unit | "concn (mmol/mol)": per mol amino acid (= per mol MGO, equimolar). The text converts 178.8 mmol/mol to "∼17.9 mol %" |
| Replicates | triplicates (RSD printed per cell) |
| Cocoa | Sulawesi (SU) and Ghana (GH) beans: unfermented sun-dried (UF), 7 days fermented (7F), each roasted (-R): 15 min at 95 °C, then raised to 115 °C within 20 min at 1 °C/min, pilot roaster. Amines by dansyl-LC-MS/MS; aldehydes by SIDA after SAFE, 2D-GC-MS; free amino acids by amino-acid analyser (norleucine IS) |

## 4. Results [M]

### 4.1 Tables 6-9 (MGO + amino acid, 170 °C; mmol/mol amino acid, RSD % in parentheses)

**Table 6, L-Leu: 3-methylbutanal / 3-methylbutylamine**

| min | aldehyde A | aldehyde B | amine A | amine B |
|---|---|---|---|---|
| 2.5 | 27.7 (10.2) | 18.0 (22.5) | 0.02 (18.1) | 0.02 (18.2) |
| 5 | 36.8 (14.1) | 28.5 (2.2) | 0.04 (2.6) | 0.03 (22.6) |
| 10 | 57.6 (12.4) | 48.7 (10.7) | 0.20 (9.2) | 0.08 (9.4) |
| 20 | 149.4 (6.5) | 64.6 (2.6) | 0.93 (8.4) | 0.12 (24.3) |
| 30 | 178.8 (4.5) | 84.6 (3.3) | 3.44 (1.8) | 0.26 (14.4) |

**Table 7, L-Ile: 2-methylbutanal / 2-methylbutylamine**

| min | aldehyde A | aldehyde B | amine A | amine B |
|---|---|---|---|---|
| 2.5 | 67.8 (15.2) | 31.8 (12.7) | 0.04 (21.6) | 0.02 (7.4) |
| 5 | 78.0 (12.2) | 41.2 (9.6) | 0.12 (31.8) | 0.03 (1.2) |
| 10 | 142.5 (14.1) | 68.5 (19.9) | 0.40 (31.3) | 0.05 (22.6) |
| 20 | 294.9 (2.4) | 117.6 (3.1) | 1.44 (21.9) | 0.15 (15.0) |
| 30 | 402.0 (4.1) | 122.1 (11.2) | 2.43 (6.6) | 0.13 (12.9) |

**Table 8, L-Phe: phenylacetaldehyde / 2-phenylethylamine**

| min | aldehyde A | aldehyde B | amine A | amine B |
|---|---|---|---|---|
| 2.5 | 19.4 (13.4) | 9.0 (4.7) | 0.28 (30.2) | 0.02 (5.9) |
| 5 | 23.2 (0.5) | 22.8 (5.4) | 0.51 (16.6) | 0.09 (2.7) |
| 10 | 28.4 (2.7) | 43.6 (1.7) | 1.95 (14.8) | 0.26 (4.8) |
| 20 | 26.4 (3.5) | 63.0 (5.9) | 6.53 (6.9) | 0.54 (7.4) |
| 30 | 26.2 (8.8) | 72.8 (3.2) | 6.13 (11.2) | 0.78 (15.4) |

**Table 9, L-Met: 3-(methylthio)propanal / 3-(methylthio)propylamine**

| min | aldehyde A | aldehyde B | amine A | amine B |
|---|---|---|---|---|
| 2.5 | 15.7 (5.4) | 9.1 (23.2) | 0.02 (27.6) | 0.01 (37.4) |
| 5 | 21.6 (13.4) | 15.3 (3.1) | 0.08 (12.3) | 0.09 (11.1) |
| 10 | 33.3 (10.7) | 19.3 (1.0) | 0.51 (20.9) | 0.19 (10.6) |
| 20 | 19.6 (8.7) | 25.3 (9.9) | 2.25 (0.4) | 0.45 (13.1) |
| 30 | 12.6 (15.4) | 22.5 (2.1) | 4.71 (18.9) | 0.72 (18.8) |

Text summaries [M/Q]: 3-methylbutanal "∼17.9 mol %" and its amine "∼0.3 mol %" after 30 min (model A); PEA "∼0.6 mol %"
after 30 min, PA : PEA "4:1 in model A"; methional peaks at 10 min in model A and is then degraded, the Met amine reaches
"∼0.5 mol %" (model A); Strecker amines overall "0.2-0.6 mol %"; "the formation of the respective Strecker aldehyde was
always clearly favored, in particular from L-leucine and L-isoleucine".

Shape notes [Z]:
- Model B aldehydes decelerate. PA increments are 13.8 / 20.8 / 19.4 / 9.8 mmol/mol per interval; 2-methylbutanal is near a
  plateau after 20 min.
- Model A is faster for Leu / Ile (x2 - 3 of B at 30 min) and stalls early for Phe (plateau ~26-28 from 10 min) and Met
  (decline after 10 min). That shows aldehyde loss in the low-water pot at 170 °C.
- Aldehydes are non-zero at 2.5 min in a closed vessel at 170 °C, so the heat-up is fast or reaction starts during heat-up.
  The heat-up is not reported.

### 4.2 PEA + MGO (text p. 1737, "data not shown")

Model A conditions: ~0.8 mol % phenylacetaldehyde, 7.1 mol % benzaldehyde (main product), plus formaldehyde proposed.

### 4.3 Cocoa (Tables 2, 4, 5; µg/kg, RSD of triplicates)

| sample | PEA SU / GH | 2-methylpropylamine SU / GH | 2- + 3-methylbutylamine SU / GH | 3-methylbutanal SU / GH | 2-methylbutanal SU / GH | phenylacetaldehyde SU / GH |
|---|---|---|---|---|---|---|
| UF | 371 / 227 | 57 / 56 | 256 / 129 | 58 / 116 | 109 / 143 | 42 / 16 |
| 7F | 723 / 1168 | 191 / 331 | 498 / 1219 | 1152 / 1636 | 1132 / 2075 | 73 / 34 |
| UF-R | 779 / 624 | 369 / 303 | 914 / 966 | 565 / 511 | 665 / 504 | 43 / 91 |
| 7F-R | 7746 / 10216 | 1988 / 2793 | 9865 / 17070 | 6841 / 8470 | 5718 / 3791 | 63 / 202 |

Free amino acids (Table 5, mg/kg, means of triplicates; UF / 7F / UF-R / 7F-R): Phe SU 66 / 1100 / 83 / 700, GH 190 / 1120 /
120 / 700; Val SU 95 / 510 / 120 / 470, GH 190 / 600 / 130 / 460; Leu SU 81 / 1100 / 97 / 800, GH 170 / 1240 / 100 / 760;
Ile SU 87 / 320 / 100 / 290, GH 140 / 390 / 100 / 280; Met SU 9 / 190 / 9 / nd, GH 15 / 220 / 15 / 87.
3-(Methylthio)propylamine was not detected in cocoa. Methylpropanal and methional were not measured in cocoa.

### 4.4 Molar ratios in roasted fermented cocoa [Z]

MW: 2-/3-methylbutanal 86.13, the amines 87.16, PA 120.15, PEA 121.18.

| ratio | SU | GH |
|---|---|---|
| (2- + 3-methylbutanal) / (2- + 3-methylbutylamine) | 1.29 | 0.73 |
| PEA / PA | ~120 | ~50 |

The authors' reading [Q]: "nearly identical concentrations" of the C5 amines and aldehydes in cocoa, and much more PEA than
PA. Fermentation is needed to make an amine precursor. Roasting did not remove Phe faster than Leu, so a lack of Phe cannot
explain the low PA. They propose "thermogenic amine" for amines formed this way.

## 5. Flags

- [!] **Model amino acids are Ile / Leu / Met / Phe**, not Val. No methylpropanal or 2-methylpropylamine model data exist.
  The repository's description (MECHANISM_REVIEWS sec. 3 table; the strecker review lists "Ile, Leu, Phe and Val") should be
  corrected.
- [!] **170 °C only, closed vessels; no heat-up, headspace or atmosphere stated.** The authors say outright that the
  conditions "were not meant to simulate the roasting conditions of cocoa".
- [!] **Aldehyde only, no Strecker acid.** In air at 98 °C, Hofmann 2000 finds PAA ~1.6 - 2.3x PA from MGO, so aldehyde
  yields here are not the total Strecker flux.
- [!] **Amines other than PEA are not quantified by SIDA** (butylamine IS with response factors). Their RSDs reach 37 %. The
  sub-0.1 mmol/mol amine values sit near their quantitation limit (not given in mmol/mol).
- [!] Two detection-limit statements: LoQ 0.2 µg/kg (p. 1733) vs "limit of detection ... 0.1 µg/kg" (p. 1734). They are
  probably different quantities; harmless.
- [!] Basis "mmol/mol" is not defined in the table headers; per amino acid is the natural reading. MGO and AA are equimolar,
  so either basis gives the same number.
- Cocoa rows: real-food concentrations with uncontrolled precursors (fermentation products, Amadori compounds). They are
  matrix comparators at most, and give no kinetics.

## 6. Supplementary material

**None.** The paper has no "Supporting Information Available" paragraph and no SI reference. Two results are "data not
shown" and are lost: the PEA + MGO product yields beyond the two percentages printed, and butylamine / alkylamines in
cocoa. Nothing to fetch.

## 7. Keying recommendation

- **Tables 6-9 are worth two-reader keying, but as HOLD-OUTS / signatures, not as step-2 fit rows.** That is 80 values with
  RSDs; the text layer is clean (verified against the raster), ~1 h for two readers.
  - Roles: (i) the amine-exit signature (amine share <= 1 % of aldehyde for Leu / Ile / Phe in water, pH 7). (ii) The donor-level
    identity ratios Ile : Leu : Phe : Met with MGO, a cross-check on f_AA.
  - (iii) The 170 °C model-B aldehyde courses, a hold-out for the temperature law (roadmap step 3), predicted from a 98 °C fit
    plus an Ea. The unknown heat-up and atmosphere set a declared floor of ~0.1 dex on those rows.
  - Model A (silica, 10 % water) rows: low-moisture diagnostic only.
- The cocoa tables are not worth keying for the kinetic lanes.
