# Snel 2023 (PhD thesis) — EXTRACTION (Processing and flavoring meat analogues: the aldehyde and ketone/ester binding chapters, their raw headspace plots, and the general discussion's extrapolation to extrusion; settles the K_ald unit at 10⁻¹ L/g)

**Source on disk:** `data/articles/Snel2023b.pdf` (the Wageningen thesis PDF, 178 PDF pages, text layer
present; equations are images without text). Read by eye on 2026-10-09: the table of contents, then only
Chapter 5 (ketones and esters, PDF p. 99-123), Chapter 6 (aldehydes, PDF p. 125-147), Chapter 7 (general
discussion, PDF p. 149-161) and the propositions (PDF p. 2). Chapters 1-4 (extrusion, rework, pectin) were
not read. Every number below was read from the page image and cross-checked with `pdftotext -layout`; the
raw-headspace figures (Fig. 5.2a, 6.2a/c/i, 7.2) were re-rendered at 220-300 dpi and read by eye. **Page
numbers below are the thesis's printed page numbers; PDF page = printed + 2** (printed 134 is PDF 136).
Written to settle the unit question left open by `snel2023_extraction.md` (the LWT aldehyde paper, which is
Chapter 6 verbatim) and to record what the thesis adds for flavouring plant-based meat.

| field | value |
|---|---|
| Title | "Processing and flavoring meat analogues" |
| Author | Silvia J. E. Snel |
| Promotors | Prof. Dr A.J. van der Goot (Wageningen University & Research); Prof. Dr M. Beyrer (HES-SO Valais) |
| Venue | PhD thesis, Wageningen University, defended 3 November 2023; "176 pages" (colophon); ISBN 978-94-6447-766-5 |
| DOI | 10.18174/633616 |
| Chapters used | Ch. 5 = Snel et al. 2023 Heliyon 9(6):e16503; Ch. 6 = Snel et al. 2023 LWT 185:115177; Ch. 7 = unpublished general discussion |

## 1. Methods

**Same rig as the papers.** Commercial isolates SPI (Supro 500E A), PPI (Nutralys F85M), FBPI
(FFBP-90-C-EU), CPPI (FCPP-70), WPI (BiPRO). Isolate dispersions 5, 10, 20, 30, 50 g isolate/kg in
demineralised water (unbuffered), converted to protein concentration c_p (g protein/kg) with the protein
content, which Table 5.2 (p. 106) prints as **"g/100 g dw"**, and the moisture content (§5.4.2, p. 108:
"3.1-4.6 g/kg to 31-46 g/kg"). In Ch. 5 the stock was stirred 1 h at 21 °C, Ultra-Turraxed and hydrated
24 h at 5 °C (§5.3.5, p. 104); all flavoured dispersions equilibrated **24 h at 21 °C**. Static headspace,
5 mL into APCI-Q-TOF (Xevo G2-XS) through a Venturi interface, independent triplicates. RHC % = peak area
in dispersion / peak area in water × 100 (Eq. 5.12, p. 104; Eq. 6.8, p. 129). Flavour concentrations:
Table 5.1 (p. 102) and Table 6.2 (p. 129), identical to the papers.

**Model (Eq. 5.14, p. 105; Eq. 6.6, 6.7 p. 127; Eq. 6.9, p. 130).** c_fg/c^p_fg = 1 + K_p·c_p with
K_p = a_p·P_ow (ketones, esters) or a_p·P_ow + K_ald (K_alk for 2-alkenals). The thesis prints no unit
inside the equations; c_p carries **g kg⁻¹** on every RHC and fit axis (Fig. 5.2, 5.3, 6.2, 6.3), and
a_p carries **"10⁻⁵ L/g"** in Table 5.3 (p. 111) and **"10⁻⁴ L/g"** in Table A.2 (p. 118-119). Chapter 6
(p. 130) quotes the same a_p as "4.8E-5, 1.1E-4, 8.6E-5, 1.7E-4, and 7.2E-5 g/L" ("g/L" sic, as in the
paper), and Eq. 6.9 carries the same typesetting slip as the paper (RHC set equal to its reciprocal). In Ch.
6 a_p is frozen at the ester values and only K_ald/K_alk is fitted (SciPy), one protein × one aldehyde.

## 2. Findings that matter

### 2.1 The unit question: the raw headspace data settle it at 10⁻¹ L/g

**What the thesis prints.** Table 6.3 (p. 134-135) has the same header as the paper, "K_ald (10⁻² L/g)"
and "K_alk (10⁻² L/g)", and the same 40 K values to the digit. The Ratio column differs from the LWT
paper's in ten cells by 1-3 units (thesis / paper: FBPI butanal 569/568, SPI butenal 1758/1755, SPI
hexenal 237/240, PPI butenal 295/294, PPI hexenal 57/58, FBPI butenal 476/475, FBPI hexenal 21/22, CPPI
butenal 1709/1706, WPI butenal 122/121, WPI hexenal 40/41); both sets reproduce only on the 10⁻¹ reading
(SPI butenal: 0.340 / (4.8e-5 × 10^0.60) = 1779, derived here). The Discussion worked numbers (p. 140:
100x decanal reduction at "48, 30, 8, 28, and 93 g/kg", "roughly 800, 1300, 5300, 1400, and 400x" at
"400 kg/kg" sic) are the paper's, and match the 10⁻¹ reading (`snel2023_extraction.md` §2.2).

**What the paper did not show: the raw RHC points against c_p.** Figure 6.2 (p. 132-133) plots the
measured K^eff/K^f (= RHC/100) against c_p in g/kg for every protein × aldehyde. The model's slope is
(100/RHC − 1)/c_p = K_p in kg/g ≈ L/g, so the raw points fix the absolute unit independently of every
printed constant. Points read from graph, approx. (c_p also read from the axis); least-squares slope
through the origin of 100/RHC − 1 against c_p, derived here:

| panel | c_p (g/kg) | RHC % read | slope (L/g) | K_p on 10⁻¹ reading | K_p on header 10⁻² reading |
|---|---|---|---|---|---|
| SPI butanal (6.2a) | 4, 8, 15, 22, 38 | 97, 86, 81, 72, 62 | 0.0165 | 0.0172 | 0.0019 |
| SPI hexanal (6.2a) | same | 90, 78, 72, 64, 39 | 0.036 | 0.0370 | 0.0064 |
| SPI octanal (6.2a) | same | 47, 36, 24, 19, 13 | 0.186 | 0.184 | 0.044 |
| PPI butanal (6.2c) | 3, 7, 14, 21, 34 | 89, 78, 73, 65, 53 | 0.0265 | 0.0264 | 0.0030 |
| PPI hexanal (6.2c) | same | 89, 86, 77, 78, 30 | 0.049 | 0.0529 | 0.0115 |
| PPI octanal (6.2c) | same | 40, 30, 20, 16, 11 | 0.250 | 0.254 | 0.085 |
| WPI butanal (6.2i) | 5, 9, 18, 28, 46 | 90, 81, 76, 71, 45 | 0.023 | 0.0263 | 0.0029 |
| WPI hexanal (6.2i) | same | 95, 83, 80, 70, 60 | 0.015 | 0.0155 | 0.0056 |
| WPI octanal (6.2i) | same | 55, 42, 29, 23, 16 | 0.119 | 0.118 | 0.051 |

K_p = a_p·10^logP + K_ald with a_p as printed in L/g, log P from Table 6.2, K_ald from Table 6.3 (derived
here). Nine of nine fits land within 15 % of the 10⁻¹ reading and 2.3-8.8x above the header reading. Read
directly: SPI octanal at 38 g/kg sits at about 13 % (graph); the 10⁻¹ reading predicts 12.5 %, the header
reading 37 % (derived here). The fitted lines in Fig. 6.3 (p. 136-137) run through these same points
(y-axis to 500), as the dossier on the paper already found for SPI decenal.

**The ketone/ester chapter passes the same test, so a_p is in L/g and c_p in g/kg.** Figure 5.2a (p. 109,
SPI ketones, read from graph, approx.): decanone 51, 35, 21, 15, 9 % and octanone 88, 76, 65, 57, 48 % at
c_p 3.8, 7.5, 15.0, 22.5, 37.5 g/kg (tick labels as printed). With a_p = 16 × 10⁻⁵ L/g (Table 5.3) the model
gives decanone 50.9, 34.5, 20.8, 14.9, 9.5 % and octanone 90.8, 83.4, 71.5, 62.6, 50.1 %; with a_p ten times
smaller, decanone 91 → 51 % (derived here). The text (p. 109) prints decanone RHC "4.5-51% to 0.7-9.1%"
across proteins; the model gives SPI 51.0 → 9.4 % and CPPI 6.5 → 0.7 % (derived here). The chemical leg
(a_p in L/g, c_p in g protein/kg) is therefore anchored to raw data in both chapters.

**The thesis contradicts itself once, in the other direction.** Table 7.1 (p. 152-153, general
discussion) predicts RHC at 32 % protein (320 g/kg). Its ketone and ester columns reproduce the L/g a_p,
and its aldehyde columns reproduce the **printed header** 10⁻² L/g, in all 16 cells (derived here):

| RHC % at 320 g/kg | hexanal | hexenal | octanal | octenal |
|---|---|---|---|---|
| SPI printed / 10⁻² / 10⁻¹ | 33 / 32.7 / 7.8 | 6 / 6.4 / 0.7 | 7 / 6.6 / 1.7 | 2 / 1.6 / 0.2 |
| PPI | 21 / 21.3 / 5.6 | 10 / 9.8 / 1.2 | 4 / 3.5 / 1.2 | 2 / 2.2 / 0.3 |
| FBPI | 11 / 10.5 / 1.4 | 23 / 23.2 / 4.0 | 2 / 1.6 / 0.2 | 6 / 5.6 / 1.3 |
| CPPI | 20 / 20.2 / 10.5 | 15 / 14.8 / 2.5 | 2 / 2.3 / 0.7 | 0 / 0.3 / 0.0 |

Ketone/ester columns, printed (model with Table 5.3 a_p in L/g, derived here): SPI hexanone 53 (52.9),
hexanoate 23 (22.9), octanone 11 (10.5), octanoate 3 (3.0); CPPI 6 (5.8), 8 (7.8), 1 (0.6), 1 (0.9).

**Verdict.** The absolute unit of Table 6.3 is **10⁻¹ L/g** (equivalently, the printed numbers × 0.1 give
L/g, c_p in g protein/kg). Evidence, in order of weight: (1) the raw RHC points of Fig. 6.2 give slopes
equal to K_p on the 10⁻¹ reading for nine protein × aldehyde panels and are 2.3-8.8x off on the header reading;
(2) the same raw-data test passes for Chapter 5, so a_p and c_p carry the units the model assumes; (3) the
Ratio column and the Discussion worked numbers (Ch. 6, p. 134-140). The header "10⁻² L/g" is a typesetting
error carried from the paper into the thesis, and the general discussion's Table 7.1 was computed from the
header and overstates the aldehyde headspace left by about 2-9x (above). Confidence about 95 %: the residual
doubt is the graph reading itself and the assumption that Fig. 6.2's c_p axis is the c_p used in the fit,
which the caption states.

### 2.2 What else the thesis adds

**Flavour added before high-moisture extrusion (§7.4.1, p. 152).** One unquantified test: "Some tests with
PPI and hexenal in a HME setting were performed ... no detectable concentration in the headspace was found,
even though according to our prediction in Table 7.1 we should have 21% left." Table 7.1's PPI hexenal cell
is 10 %; 21 % is the PPI hexanal cell, so either the compound or the cell is misnamed. No dose, barrel
temperature, die, sampling time or detection limit is given. The author concludes that "other factors
during processing reduce this concentration even further". On the corrected unit the prediction at 320
g/kg is 1.2 % (hexenal) or 5.6 % (hexanal) (derived here), close to what a 1-10 mg/kg dose would leave near
the APCI detection limit, so the test does not by itself show an extra processing loss. Proposition 2
(PDF p. 2): "The addition of flavors to plant proteins while extruding is not efficient."

**Heat, discussed but not measured (p. 152-153).** Extrusion at "around 140 °C"; the author expects heating
to raise covalent aldehyde binding, especially for small aldehydes, and to raise then lower hydrophobic
binding (unfolding, then aggregation), citing Wang & Arntfield 2015 (hexanal on canola, 95 °C, 60 min) and
Kühn et al. 2008 (nonanal on whey, 80 °C). No temperature other than 21 °C was measured for binding.

**Soluble versus insoluble protein, buffered (§7.4.2, Fig. 7.2, p. 154).** New data: protein fractions
separated by two-step centrifugation (16,000 × g, 30 min) in 0.01 M phosphate **pH 8**, at 0.5 and 2 %
protein; y-axis "Retention (%)" (the caption says relative headspace concentration). Read from graph,
approx.:

| retention % | soluble 0.5 / 2 % | insoluble 0.5 / 2 % | total 0.5 / 2 % |
|---|---|---|---|
| PPI methyl hexanoate | 5 / 8 | 8 / 35 | 3 / 32 |
| PPI methyl octanoate | 36 / 66 | 64 / 91 | 57 / 92 |
| PPI 2-hexanone | ~0 / 3 | 2 / 2 | ~0 / 1 |
| PPI 2-octanone | 14 / 34 | 29 / 23 | 24 / 49 (error bars span most of the axis) |
| SPI methyl hexanoate | 1 / 9 | 13 / 34 | 8 / 30 |
| SPI methyl octanoate | 27 / 53 | 54 / 81 | 43 / 77 |
| SPI 2-hexanone | ~0 / ~0 | ~0 / ~0 | ~0 / ~0 |
| SPI 2-octanone | 28 / 12 | 21 / 7 | 13 / 13 (very wide error bars) |

The insoluble fraction binds the esters more than the soluble one. Cross-check against Ch. 5 (derived here,
assuming "%" = g protein/100 g): methyl octanoate total retention predicted 53 / 82 % (PPI) and 33 / 67 %
(SPI) at 5 / 20 g/kg, against 57 / 92 and 43 / 77 read: the unbuffered-water a_p carries over to pH 8
within about 10 points.

**Aldehyde does not change texture (§7.4.3, Fig. 7.3, p. 155).** 0.1 % butanal in 45 % dry-matter PPI,
strain sweep at 30 °C and after 140 °C → 30 °C: G′ and G″ unchanged (G′ about 2 × 10² kPa at 30 °C and
3 × 10² kPa after heating, G″ about 5-6 × 10¹ kPa; read from graph, approx.).

**Ketone and ester constants (Table 5.3, p. 111; Table A.2, p. 118-119).** a_p in 10⁻⁵ L/g, ketones / esters
(± printed): soy 16 ± 0.1 / 4.8 ± 0.2; yellow pea 25 ± 0.1 / 11 ± 0.3; faba 23 ± 0.2 / 8.6 ± 0.2; chickpea
290 ± 5.7 / 17 ± 0.4; whey 22 ± 0.1 / 7.2 ± 0.1. The per-compound fits of Table A.2 (10⁻⁴ L/g) agree with
the pooled values (e.g. methyl octanoate: SPI 0.48, PPI 1.12, FBPI 0.87, CKPI 1.68, WPI 0.71; decanone:
SPI 1.64, CKPI 29.16), so the two tables use one unit consistently; Table A.2's "Uncertainty a_p" column
(e.g. 2.28E+01 against a value of 0.00) is not interpretable as printed. Butanone and methyl butanoate
show no retention (slope about 0; RHC rises with protein, read as a "pushing out" effect, p. 108-109);
methyl decanoate was excluded as non-linear. The text says hexanone RHC "hardly decreased (from 96-100% to
93-100%)" (p. 109), but Fig. 5.2a shows SPI hexanone falling to about 76 % (read from graph, approx.).
The ketone/ester chapter is consistent with the aldehyde chapter in a_p and c_p (same values, same units);
the inconsistency is confined to the K_ald/K_alk header.

**Other.** Sulfur: no thiol or sulfur volatile was measured anywhere; the only sulfur content is the
correlation of K_alk with Met and Cys contents (Table A.1, p. 142-143; five proteins). pH: only the isolate
pH (Table 5.2, Table 6.1) and the pH 8 buffer of Fig. 7.2. Cross-reference slips: §7.4.2 cites "Figure
7.1" for Fig. 7.2, §7.4.3 "Figure 7.2" for Fig. 7.3, and §7.5 puts pectin in "Chapter 6" (it is Ch. 4).

### 2.3 The author's view of what limits flavouring (§7.4, 7.5, 7.6, p. 151-157)

- Protein concentration dominates: model mixtures were 0.5-5 % isolate against about 32 % protein in a
  meat analogue (40 % isolate DM at 80 % protein). "To achieve the desired flavor concentration in a CPPI
  meat analogue, a 100x higher concentration should be added for octanone and octanoate, which is probably
  not feasible" (p. 152).
- The solid, heated matrix is untested; other ingredients may change binding via pH or salt (p. 151).
- The rotating die does not solve binding; retention is expected in any structuring process (p. 157).
- Remedies she proposes: add flavour after structuring (vacuum marination, which "mainly works for small
  meat analogue pieces"), cleaner-tasting proteins that need no maskers (breeding, extraction), and "the
  development of precurses [sic] that form flavor at the process conditions used to make the fibrous
  products" (p. 157).

## 3. What it means for the model

**The binding form is the engine's.** The matrix layer's per-gram constant K_g = (K_water/K_matrix − 1) /
protein_g_per_L (`src/kinetic_core/parameters_matrix.py`, line 313 comment) is Snel's K_p, and Fig. 6.2's
raw points measure it directly. With the unit settled, the five Snel hexanal rows enter as total K_p on
the 10⁻¹ reading, in L/g: SPI 0.0370, PPI 0.0529, FBPI 0.217, CPPI 0.0267, WPI 0.0155 (from
`snel2023_extraction.md` §2.3; spot-checked here for SPI and PPI). Live rows they would join:
`kg_hexanal_dairy` = 1.151e-2 L/g (Meynier 2002, skim milk, 30 °C) and `kg_hexanal_pea` = 2.537e-1 L/g
(Bi 2022, pea isolate, 37 °C, pH 7.6). The bracketed header-reading values in the paper's dossier can be
dropped; the class-pooling move it computed for the 10⁻¹ reading (0.054 → 0.047 L/g) is the relevant one,
and still belongs to a pre-registered wave because it moves FIT rows. The engine's matrix layer keeps
its rule "No log P term of any kind" (`parameters_matrix.py` line 50), so only total K_p enters, not the
a_p·P_ow / K_ald split; the alkenal rows stay out under the existing Michael-acceptor quarantine.

**The chain-length slope.** Live `CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81` (`parameters_matrix.py` line 626,
used at `matrix_oav.py` line 741 to derive the branched-alkanal surrogate). The unit verdict does not touch
it (a slope within one protein is unit-free); Snel's C6-C10 slopes of 2.73-3.39 per CH₂ on the 10⁻¹
reading stand as corroboration, and the C4-C6 flattening (0.77-1.47 per CH₂) stands as the warning on the
Strecker-aldehyde surrogate (paper dossier §3).

**For the purpose (meat aroma in plant products).** Use the corrected numbers, not the thesis's Table 7.1,
for any extrapolation to analogue protein loads: at 320 g protein/kg the equilibrium headspace left for
C6-C8 aldehydes is about 0-11 % of the water value on the four plant isolates, not 0.3-33 %. A Maillard
model that generates odorants in situ (the author's "precursors that form flavour at the process
conditions") faces the same partition at the end: whatever forms in the matrix is bound by the same K_p.
The thesis gives nothing on how that partition shifts at extrusion temperature.

## What it does not give

- Any binding constant at a temperature other than 21 °C, or in a heated, extruded or solid matrix.
- Retention or loss numbers for flavour added before or during extrusion: one undetected PPI-hexenal (or
  hexanal) HME headspace, with no conditions, dose or detection limit.
- Any thiol, sulfide, disulfide, pyrazine or other Maillard odorant; sulfur enters only as Met/Cys
  correlations over five proteins.
- Tabulated raw RHC data: the points exist only as figures (values above are graph readings).
- Binding at any pH other than the isolates' own (unbuffered) and the pH 8 fractionation experiment.
- An independent a_p for aldehydes (frozen from esters), or any test of reversibility for the "covalent"
  term.
- Sensory data or odour thresholds.
