# Snel, Pascu, Bodnár, Avison, van der Goot & Beyrer 2023 — EXTRACTION (C4-C10 n-alkanals and trans-2-alkenals against soy, pea, faba, chickpea and whey isolates: static headspace APCI-TOF-MS, 21 °C, Harrison & Hills partition model with a fitted covalent term)

**Source on disk:** `data/articles/Snel2023.pdf` (the publisher's PDF, 8 pages, open access CC BY-NC-ND, with a
text layer). Read and checked by eye on 2026-10-09: every number below was read from the page image and
cross-checked against `pdftotext -layout`; Table 1 (p. 3) and Table 3 (p. 6) were also re-rendered at 220 dpi
and re-read. Figures 2 and 3 (p. 4-5) have no data table; the few values taken from them are labelled "read
from graph, approx.". The supplementary material (background spectra, amino-acid correlation table) is not
on disk. Written for the matrix/OAV layer's aldehyde binding (plant-protein binding of the lipid-oxidation
off-notes and of the aldehydes a meat flavour carries), not for the core fit.

**Companion paper.** This is the aldehyde sequel of Snel et al. 2023, Heliyon 9(6):e16503 (ketones and
esters, same five isolates, same rig). That paper has no dossier in this directory; its ketone and ester
a_p values are carried in `data/lit/binding_constants.yml` (source `snel_2023_heliyon`, records
`snel2023_ap_*`, verified from the Europe PMC full text on 2026-08-27). That file is a literature ledger:
`src/data_paths.py` names it (`BINDING_CONSTANTS`) and no engine code reads it.

| field | value |
|---|---|
| Title | "Flavor-protein interactions for four plant protein isolates and whey protein isolate with aldehydes" |
| Authors | Silvia J.E. Snel, Mirela Pascu, Igor Bodnár, Shane Avison, Atze Jan van der Goot, Michael Beyrer (HES-SO Valais, Sion; Wageningen University; Firmenich, Geneva) |
| Venue | LWT - Food Science and Technology 2023, 185, 115177; received 21 Apr 2023, accepted 8 Aug 2023 |
| DOI | 10.1016/j.lwt.2023.115177 |

## 1. Methods

**Proteins (§3.1, p. 2-3; Table 1, p. 3).** SPI Supro 500E A (Solae); PPI Nutralys F85M (Roquette, the
same product as the Heliyon ketone records); FBPI FFBP-90-C-EU and CPPI FCPP-70 (AGT Foods); WPI BiPRO
(Davisco). Table 1, means ± SD (amino acids in duplicate, the rest in triplicate):

| | Soy | Yellow pea | Faba bean | Chickpea | Whey |
|---|---|---|---|---|---|
| WHI | 147.35 ± 0.63 | 149.41 ± 0.09 | 143.22 ± 1.69 | 152.28 ± 0.91 | 236.18 (from Amagliani 2017) |
| Protein (g/100 g) | 83.25 ± 2.77 | 77.11 ± 0.29 | 80.42 ± 2.99 | 67.03 ± 1.83 | 89.7 |
| pH | 7.1 ± 0.01 | 7.5 ± 0.01 | 6.4 ± 0.02 | 6.6 ± 0.01 | 7.1 ± 0.02 |
| Solubility (g/g) | 0.59 ± 0.01 | 0.41 ± 0.03 | 0.12 ± 0.01 | 0.12 ± 0.00 | 1.02 ± 0.00 |
| Moisture (g/100 g) | 8.8 ± 0.0 | 8.1 ± 0.0 | 8.0 ± 0.1 | 6.6 ± 0.1 | 5.2 ± 0.1 |

The pH is listed as a property of the isolate; how it was measured (dispersion strength) is not printed,
and the pH of the flavoured dispersions themselves is not reported. The dispersions are unbuffered
(demineralised water).

**Dispersions (§3.2, p. 3).** Stock 50 g/kg, diluted to 5, 10, 20, 30, 50 g isolate/kg in demineralised
water, then converted to protein concentration c_p (g protein/kg) with Table 1's protein and dry matter.
Each aldehyde pre-diluted in ethanol (final ethanol 1 g/kg, "did not affect the MS signal"), added,
vortexed 30 s, **equilibrated 24 h at 21 °C**. Flavour concentrations (Table 2, p. 3, mg/kg): butanal 1.00,
hexanal 0.75, octanal 0.25, decanal 0.05, trans-2-butenal 5.00, trans-2-hexenal 1.00, trans-2-octenal
0.75, trans-2-decenal 0.50. Table 2's log P (EPIWEB v4.1 KOWWIN v1.68): 0.60, 1.80, 2.78, 3.76 (alkanals
C4-C10) and 0.60, 1.58, 2.57, 3.55 (alkenals). Two printed inconsistencies: §3.2 says the range is "0.05 to
1.00 mg/kg" while Table 2 gives butenal 5.00; §4.1 (p. 4) says "0.05 to 1.00 g/kg".

**Measurement (§3.3, p. 3).** Static headspace, 5 mL injected into an APCI G2-XS Q-TOF (Waters) through a
Venturi interface, independent triplicates. Relative headspace concentration RHC % = (peak area flavoured
dispersion − peak area blank dispersion) / peak area flavour in water × 100 (Eq. 8). No internal standard;
the water leg is the reference, so every number is a within-run ratio.

**Model (§2, §3.4, p. 2-3).** Harrison & Hills / Viry: K_p = a_p·P_ow + K_ald (Eq. 7; K_alk for alkenals),
fitted through Eq. 9, c_fg/c^p_fg = 1 + (a_p·P_ow + K_ald)·c_p (Fig. 3 plots this, water headspace over
dispersion headspace, against c_p in g/kg). As printed, Eq. 9 sets this equal to RHC = K^eff/K^f, which is
its reciprocal; the fits in Fig. 3 are of the reciprocal, so this is a typesetting slip. a_p is NOT fitted:
it is fixed at the ester values of the Heliyon paper, printed as "4.8E-5, 1.1E-4, 8.6E-5, 1.7E-4, and
7.2E-5 g/L for SPI, PPI, FBPI, CPPI, and WPI" (p. 3; "g/L" sic, the Heliyon unit is L/g; the soy and pea
values match `binding_constants.yml`'s Heliyon ester records, 4.8 and 11 × 10⁻⁵ L/g). Only K_ald / K_alk
is fitted (SciPy), one protein × one aldehyde at a time. Ratio = K_ald / (a_p·P_ow) (Eq. 10). Statistics:
ANOVA and Tukey in R; R² reported as the squared Pearson correlation.

## 2. Findings that matter

### 2.1 Table 3 (p. 6), as printed: K_ald and K_alk, header unit "10⁻² L/g", fit ± uncertainty (n = 3)

| protein | aldehyde | K_ald | R² | Ratio | alkenal | K_alk | R² | Ratio |
|---|---|---|---|---|---|---|---|---|
| SPI | butanal | 0.17 ± 0.01 | 0.78 | 87 | butenal | 3.40 ± 0.30 | 0.88 | 1755 |
| | hexanal | 0.34 ± 0.03 | 0.84 | 11 | hexenal | 4.40 ± 0.23 | 0.91 | 240 |
| | octanal | 1.55 ± 0.05 | 1.00 | 5 | octenal | 17.13 ± 0.25 | 0.99 | 96 |
| | decanal | 17.88 ± 1.29 | 0.91 | 6 | decenal | 95.64 ± 2.87 | 0.98 | 56 |
| PPI | butanal | 0.26 ± 0.02 | 0.86 | 58 | butenal | 1.32 ± 0.05 | 0.99 | 294 |
| | hexanal | 0.46 ± 0.08 | 0.73 | 7 | hexenal | 2.46 ± 0.05 | 0.99 | 58 |
| | octanal | 1.88 ± 0.06 | 1.00 | 3 | octenal | 9.81 ± 0.12 | 0.99 | 24 |
| | decanal | 27.03 ± 2.28 | 0.93 | 4 | decenal | 60.95 ± 1.36 | 0.98 | 15 |
| FBPI | butanal | 1.96 ± 0.13 | 0.94 | 568 | butenal | 1.64 ± 0.17 | 0.78 | 475 |
| | hexanal | 2.12 ± 0.16 | 0.85 | 39 | hexenal | 0.71 ± 0.07 | 0.69 | 22 |
| | octanal | 13.56 ± 1.56 | 0.55 | 26 | octenal | 2.09 ± 0.05 | 0.98 | 7 |
| | decanal | 126.94 ± 13.68 | 0.60 | 25 | decenal | 27.40 ± 1.70 | 0.96 | 9 |
| CPPI | butanal | 0.13 ± 0.03 | 0.10 | 20 | butenal | 11.46 ± 1.39 | 0.76 | 1706 |
| | hexanal | 0.16 ± 0.02 | 0.96 | 1 | hexenal | 1.15 ± 0.03 | 0.98 | 18 |
| | octanal | 3.24 ± 0.10 | 0.98 | 3 | octenal | 112.86 ± 58.56 | 0.01 | 183 |
| | decanal | 25.45 ± 1.04 | 0.98 | 3 | decenal | 429.36 ± 62.14 | 0.64 | 72 |
| WPI | butanal | 0.26 ± 0.04 | 0.55 | 91 | butenal | 0.35 ± 0.02 | 0.84 | 121 |
| | hexanal | 0.11 ± 0.01 | 0.88 | 2 | hexenal | 1.10 ± 0.04 | 0.96 | 41 |
| | octanal | 0.75 ± 0.03 | 1.00 | 2 | octenal | 3.21 ± 0.06 | 0.99 | 12 |
| | decanal | 6.51 ± 1.01 | 0.68 | 2 | decenal | 14.56 ± 0.15 | 1.00 | 6 |

Fits to flag: CPPI butanal (R² 0.10), CPPI octenal (R² 0.01, ± 52 %), CPPI decenal (0.64), FBPI octanal and
decanal (0.55, 0.60), WPI butanal and decanal (0.55, 0.68). The authors attribute the CPPI alkenal
uncertainty to headspace near the APCI detection limit (§4.2, p. 6).

### 2.2 The header unit is contradicted by the paper's own arithmetic, by exactly 10x (derived here)

Recomputing the Ratio column from Table 3, the printed a_p and Table 2's log P with K in 10⁻² L/g gives one
tenth of every printed ratio, in all 40 cells: SPI butanal 8.9 vs 87 printed (0.0017 / (4.8e-5 × 10^0.60
= 1.91e-4)), SPI hexanal 1.1 vs 11, PPI decanal 0.43 vs 4, SPI butenal 178 vs 1755, CPPI hexanal 0.15 vs 1.
With K in 10⁻¹ L/g every printed ratio is reproduced to rounding.

The Discussion's worked examples (§5, p. 7) say the same thing independently. "A 100x reduction of the
headspace concentration of decanal is already reached at a protein concentration of 48, 30, 8, 28, and
93 g/kg for SPI, PPI, FBPI, CPPI, and WPI"; for decenal "10, 15, 33, 2, and 58 g/kg". Solving
1 + (a_p·P_ow + K)·c_p = 100:

| | decanal, K in 10⁻² L/g | decanal, K in 10⁻¹ L/g | printed | decenal, 10⁻² | decenal, 10⁻¹ | printed |
|---|---|---|---|---|---|---|
| SPI | 218 | 48.0 | 48 | 88 | 10.2 | 10 |
| PPI | 110 | 29.7 | 30 | 99 | 15.3 | 15 |
| FBPI | 56 | 7.5 | 8 | 171 | 32.5 | 33 |
| CPPI | 80 | 28.1 | 28 | 20 | 2.3 | 2 |
| WPI | 207 | 92.9 | 93 | 247 | 57.8 | 58 |

The "400 kg/kg protein" examples (sic; 400 g/kg is meant) match the same reading: decanal reductions
827, 1335, 5277, 1410, 427x (printed "roughly 800, 1300, 5300, 1400, and 400x"), decenal 3895, 2595, 1219,
17417, 686x (printed 3900, 2600, 1200, 17400, 700x). Fig. 3 agrees too: the SPI decenal line reaches
about 490 at c_p = 50 g/kg (read from graph, approx.; 1 + 9.73 × 50 = 488 on the 10⁻¹ reading, 57 on the
header reading), and PPI octanal's RHC is "around 40-50%" at 3-4 g/kg in §4.1 (10⁻¹ reading at 3.5 g/kg:
53 %; header reading: 77 %).

**Conclusion.** Three independent internal checks (ratio column, Discussion examples, Fig. 3 fit lines)
agree that the Table 3 numbers are in units of 10⁻¹ L/g (with a_p in L/g and c_p in g protein/kg), and
none supports the printed 10⁻² L/g. Either the header is wrong, or the header is right and the ratios,
examples and figures were all computed with K ten times larger. The paper cannot be used for absolute
constants until the authors (or the "on request" data) settle it. Every derived number below is given
on the 10⁻¹ reading, with the header reading in brackets.

### 2.3 Total per-gram constants, K_p = a_p·P_ow + K (derived here, L/g)

| | butanal | hexanal | octanal | decanal | hexenal | octenal | decenal |
|---|---|---|---|---|---|---|---|
| SPI | 0.0172 [0.0019] | 0.0370 [0.0064] | 0.184 [0.044] | 2.06 [0.455] | 0.442 [0.046] | 1.73 [0.189] | 9.73 [1.13] |
| PPI | 0.0264 [0.0030] | 0.0529 [0.0115] | 0.254 [0.085] | 3.34 [0.903] | 0.250 [0.029] | 1.02 [0.139] | 6.49 [1.00] |
| FBPI | 0.196 [0.020] | 0.217 [0.027] | 1.41 [0.187] | 13.2 [1.76] | 0.074 [0.010] | 0.241 [0.053] | 3.05 [0.579] |
| CPPI | 0.0137 [0.0020] | 0.0267 [0.0123] | 0.426 [0.135] | 3.52 [1.23] | 0.122 [0.018] | 11.4 [1.19] | 43.5 [4.90] |
| WPI | 0.0263 [0.0029] | 0.0155 [0.0056] | 0.118 [0.051] | 1.07 [0.479] | 0.113 [0.014] | 0.348 [0.059] | 1.71 [0.401] |

Chain-length slope per CH₂ of K_p, C6 to C10 (4 CH₂): SPI 2.73 [2.90], PPI 2.82 [2.97], FBPI 2.79 [2.85],
CPPI 3.39 [3.16], WPI 2.88 [3.04]. From C4 to C6 it is 0.77-1.47 per CH₂: butanal is bound almost as much
as hexanal, which the authors put down to covalent chemistry dominating at C4 (Ratio 20-568).

Same-carbon unsaturation contrast, K_p(hexenal)/K_p(hexanal): SPI 11.9 [7.1], PPI 4.7 [2.5], FBPI 0.34
[0.39], CPPI 4.5 [1.5], WPI 7.3 [2.4]. Faba bean binds hexanal more than hexenal at every chain length.

### 2.4 Other printed results

- §5, p. 7: octanal 70 % bound in a 7 g/kg PPI dispersion (from headspace), against Wang & Arntfield
  2014's 68 % at 10 g/kg.
- §4.3, p. 6: K_alk correlates with methionine (75, 43, 73, 75 % for C4-C10) and cysteine (30, 84, 26,
  29 %); K_ald correlates negatively with Cys and Met (about −80 %) and with histidine 55-59 %, arginine
  61-73 %, leucine 86-89 %, valine 51-78 %. Five proteins per correlation; the authors call it "not
  conclusive".
- Ranking (§6): aldehydes most retained by FBPI and PPI, alkenals by SPI and CPPI; WPI lowest except
  butanal.

## 3. What it means for the model

**The form is the engine's form.** The matrix layer's per-gram constant is K_g = (K_water/K_matrix − 1)/c_p
(`src/kinetic_core/parameters_matrix.py`, REVERSIBLE_BINDING comment block), which is exactly Snel's
c_fg/c^p_fg − 1 = K_p·c_p. So each Table 3 cell enters as its TOTAL K_p (section 2.3), one MatrixParameter
per protein × compound, method `static_headspace_partition`, temperature 21 °C, with a new MATRIX_LOADING
per isolate (Table 1's protein fraction, pH as printed for the isolate). Two things stay out:

- **The a_p·P_ow / K_ald split.** It is an attribution by construction (a_p frozen from esters), not a
  measurement, and the layer forbids a log P term ("No log P term of any kind", k4b hold-out guard #4,
  module docstring). The covalent label is not tested either: no reversibility arm was run, while the
  engine's `COVALENT_CEILING` (source anchor meynier2004 / shepelev2024, `parameters_matrix.py`) carries
  hexanal as `reversible_share_headspace_timescale = 0.98`. Snel's "covalent" K_ald cannot feed the
  covalent term, which in any case contributes 0.0 to every point prediction.
- **The alkenal rows**, on the precedent already applied to `kg_t_2_hexenal_dairy` and `kg_t_2_octenal_pea`
  (quarantined as binding constants: Michael acceptors on a lysine-rich protein, measured by
  disappearance). Snel's 24 h equilibration is longer than Bi's 2 h, so the quarantine applies with more
  force.

**Comparison with live values.** `fit_class_binding_constants()` (`src/kinetic_core/matrix_oav.py`, run
2026-10-09) pools `n_alkanal` at **0.0540 L/g**, the geometric mean of `kg_hexanal_dairy` 1.151e-2 L/g
(Meynier 2002, skim milk, 30 °C) and `kg_hexanal_pea` 2.537e-1 L/g (Bi 2022, pea isolate, 37 °C, pH 7.6).
Snel's PPI hexanal is 0.0529 [0.0115] L/g: 4.8x [22x] below Bi's pea row (Bi's isolate was made in the
laboratory, Snel's is the commercial Nutralys F85M; Bi ran 2 h at 37 °C and pH 7.6), and on the 10⁻¹
reading it lands on the pooled class value. Snel's WPI hexanal, 0.0155 [0.0056], sits 1.35x above [2.1x
below] Meynier's skim-milk row. Adding the five Snel hexanal rows to the two FIT rows would move the
pooled class from 0.054 to 0.047 L/g on the 10⁻¹ reading and to 0.017 L/g on the header reading
(geometric means, derived here): a small move or a 3.2x one, in a FIT-row quantity; that moves scored predictions and belongs
to a pre-registered wave, after the unit is settled.

**The chain-length slope is corroborated.** The live `CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81`
(`parameters_matrix.py`; Andriot 2.72, Damodaran 2.90, Guo 2.86) against Snel's C6-C10 slopes of
2.73-3.39 on five more proteins (both unit readings give 2.85-3.39; the slope is unit-free within a
protein). Below C6 the slope breaks down (0.77-1.47 per CH₂ from C4 to C6). The engine's
`branched_alkanal` surrogate divides the C6 constant by 2.81 to get C5 (2- and 3-methylbutanal); Snel's
C4-C6 data suggest a C5 alkanal is about 1.0-1.5x below hexanal, not 2.8x, so the surrogate likely
under-binds the Strecker aldehydes by roughly 2x (derived here; linear n-alkanals, not branched).

**The unsaturation penalty.** The live `fit_unsaturation_penalty()` returns 3.73x from two FIT rows
(Vega gelatin 2.81x, Meynier skim milk 4.95x). Snel offers five SAME-CARBON C6 pairs (the property whose
absence excluded Bi's C8/C6 pair), but they are ratios of per-gram constants, not of matrix shifts, and
they span 0.34-11.9 [0.39-7.1] across proteins measured in one lab: the penalty is protein-dependent, and
faba bean inverts it, which runs against ordinal gate G-1 (enals ≫ hexanal) for that protein.

**For a formulator.** At meat-analogue protein loads binding is not a correction, it is the whole story for
long aldehydes: on the paper's own model decanal headspace drops 400-5300x at 400 g/kg protein (§5).
Hexanal, the beany off-note, is the least bound of the series on every plant isolate (K_p 0.03-0.2 L/g),
so protein binds the off-note far less than it binds a C8-C10 fatty-meaty aldehyde added to cover it;
faba isolate is the strongest aldehyde binder (about 10x the others for K_ald), chickpea and soy the
strongest alkenal binders.

## What it does not give

- An unambiguous absolute unit for Table 3 (section 2.2). Raw RHC data are "available on request".
- Any temperature other than 21 °C; any heat-treated or extruded protein; any measured pH of the
  dispersions; any buffer.
- Any thiol, disulfide, pyrazine or other Maillard odorant: aldehydes only.
- A test of reversibility or adduct formation: "covalent" is the authors' label for what the ester-derived
  hydrophobic term does not explain.
- A free a_p for aldehydes: it is fixed from esters, so the split between the two terms carries no
  independent information.
- Odour thresholds or sensory data.
