# Wang & Arntfield 2014 — EXTRACTION (C6-C8 aldehydes and 2-ketones on salt- and alkali-extracted pea and canola isolates and wheat gluten, 1 % w/v, pH 8, 30 C / 3 h, ITEX headspace GC/MS)
### The single-flavour baseline of the Manitoba pea-binding series; every percent-bound value is in Fig. 1 and none is printed.

**Source on disk:** `data/articles/wang2014.pdf` (owner's download, 2026-09-08). Read from the
scratchpad text layer (`wang2014.txt`, clean); the PDF was re-scanned with pypdf for the word
"Table" and has none — the paper carries six figures and no table. Repo status before this dossier:
`data/lit/binding_constants.yml` already records, under source `wang_2015_umanitoba_thesis`, that
"the hexanal / heptanal / octanal x pea percent-bound values this campaign most wanted are in Figure
3.1a of this thesis and are FIGURE-ONLY"; thesis Chapter 3 is this paper, and this dossier confirms
the same verdict from the journal version. Barallat-Perez 2023 (already in `binding_constants.yml`)
describes its own headspace step as "a modified version of the Wang and Arntfield protocol" — this
is the protocol.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Binding of carbonyl flavours to canola, pea and wheat proteins using GC/MS approach" |
| Authors | Kun Wang, Susan D. Arntfield (Department of Food Science, University of Manitoba, Winnipeg) |
| Venue | Food Chemistry 157 (2014) 364-372; received 26 Nov 2013, accepted 11 Feb 2014 |
| DOI | 10.1016/j.foodchem.2014.02.042 (PII S0308-8146(14)00217-9) |
| Naming | CPIs / PPIs = salt-extracted canola / pea protein isolates (protein micellar mass route); CPIa / PPIa = alkaline-extracted, acid-precipitated. "1-heptanal" in the Fig. 3 caption = heptanal. |
| Companions | Wang & Arntfield 2015, Food Hydrocolloids 43:410-417 (mixtures + heat; dossier `wang2015_extraction.md`); Wang & Arntfield 2015, Food Res. Int. 77:1-9 (salts and pH, ketones only; the source of the `wang2015_ppi_*_ph7` records already in `binding_constants.yml`). |
| Compound registry | hexanal -> `hexanal`; heptanal -> `heptanal`; octanal -> not in registry; 2-hexanone, 2-heptanone, 2-octanone -> not in registry; 2-butyl-2-octenal, 2-pentyl-2-nonenal -> not in registry. |

## 1. Why it matters

Need (a). The engine's matrix layer (`src/kinetic_core/matrix_sites.py`) binds hexanal-class
aldehydes to a protein amine pool with rate brackets from the adduct dossiers, and
`data/species/protein_matrices.yml` lists only beta-lactoglobulin because "no dossier on disk gives
the free thiol and disulfide content of pea or soy isolates". This paper is the primary pea-isolate
measurement that the whole later literature (Barallat-Perez 2023, Bi 2022, Snel 2023) cites for
method; it measures headspace depletion of hexanal, heptanal and octanal by a lab-made pea isolate at
a stated loading (10 g powder/L, 82.68 % protein), pH 8, 30 C, after 3 h. If its numbers were printed
they would be the first `percent_bound_at_conditions` hexanal rows for pea. They are not printed
(Fig. 1a/b only), so what the repository can take is: the protocol and its conditions in full, the
protein purities, the directional orderings stated in the text, the 2-5x aldehyde-over-ketone factor,
the observation that aldehydes (not ketones) form aldol by-products on canola protein (Schiff-base
chemistry, i.e. the amine route the engine models), and the cited Lys / Cys / Met content of the pea
isolate, which is the only route on disk to an amine site density for pea (secondary source, see
§4). Need (b) is not touched: no LOX, no processing volatiles.

## 2. Methods as they matter to a model

- **Materials.** Commercial yellow pea flour (Best Cooking Pulses, Portage la Prairie); canola meal
  (BMW Canola AL018); commercial vital wheat gluten (Arrowhead Mills). Flavours analytical grade,
  Sigma-Aldrich.
- **PPIs (salt-extracted pea isolate).** Sieved flour (500 um) in 0.3 M NaCl, flour:solution 3:10
  w/v, stirred 30 min; centrifuged 4260 g / 4 C / 15 min; supernatant diluted with 2 volumes cold
  water, 3 C for 2 h (micelle precipitation); pellet at 680 g; dialysed 12-14 kDa against 20 volumes
  water, 72 h; freeze-dried.
- **PPIa (alkaline pea isolate).** 40 g flour in 600 mL water at pH 9.5 (1 M NaOH), 1 h; centrifuged
  4500 g / 20 min; supernatant to pH 4.5 (0.1 M HCl); pellet washed with pH-4.5 water; freeze-dried.
  (Canola: CPIs by PMM from 0.5 M NaCl with 10 kDa ultrafiltration; CPIa at pH 8 extraction / pH 4
  precipitation.)
- **Protein contents (Dumas, N x 5.7):** CPIs 87.32 %, CPIa 75.35 %, **PPIs 82.68 %**, **PPIa
  82.82 %**, gluten 76.01 %. So 1 % w/v powder = 10 g powder/L = **8.27 g protein/L (PPIs)** or
  8.28 g protein/L (PPIa).
- **Buffer.** 0.01 M potassium phosphate, pH 8; "The ionic strength was kept as low as possible to
  minimise the effect of salt on protein conformation." 2 % w/v protein stock dispersed by ultrasonic
  bath 20 min.
- **Flavour stock.** "Stock solutions of each volatile flavour compound were prepared in phosphate
  buffer solution at 1000 ppm (0.1 mL/100 mL)". The stock is therefore **volumetric** (1000 uL/L),
  and the working "250 ppm" is 250 uL/L. Converted with liquid densities (hexanal 0.815, heptanal
  0.818, octanal 0.821, 2-hexanone 0.812, 2-heptanone 0.820, 2-octanone 0.820 g/mL, handbook values
  not printed in the paper): **250 ppm v/v = 203-205 mg/L**, i.e. hexanal 2.03 mM, heptanal 1.79 mM,
  octanal 1.60 mM, 2-hexanone 2.03 mM, 2-heptanone 1.80 mM, 2-octanone 1.60 mM. Note that
  `binding_constants.yml` stores the sibling thesis records as `flavor_concentration_mg_per_L:
  250.0`; if the thesis stocks were made the same way that field is ~18 % high (flag 5).
- **Vial charge.** 1 mL of 2 % protein + 0.5 mL buffer + 0.5 mL flavour stock in a 20-mL crimp vial
  (2 mL liquid, 18 mL headspace); flavour added last. CPIa, PPIa and gluten were "not highly soluble
  in the desired buffer", so 0.02 g powder was weighed into the vial + 1.5 mL buffer + 0.5 mL stock
  (also 1 % w/v, but a suspension). Final: **1 % w/v isolate, 250 ppm flavour, pH 8**.
- **Equilibration.** Shaking water bath **30 C, 125 rpm, 3 h**; "Preliminary testing found that
  3 h was adequate to reach equilibrium." Duplicate vials, each sampled once.
- **Headspace sampling.** "After mixing, samples were incubated and shaken for 14 min at 40 C and
  1 mL of sample headspace was aspirated into the GC injector port by a CombiPal autosampler unit
  with PAL Itex-2 (In-Tube-Extraction) absorber attachment (CTC Analytics AG, Switzerland) after one
  absorption cycle." So the measured headspace is at **40 C**, after a 3-h 30 C equilibration: a
  dynamic (ITEX) headspace, not static, no internal standard, no calibration curve.
- **GC/MS.** Varian CP-3800 / 320-MS triple quadrupole in single-quad mode, splitless; VF-5ms
  30 m x 0.2 mm, film printed as "20 um" (almost certainly 0.2 um; flag 6); He 4 mL/min; oven ramp
  25 C/min to 265 C, hold 3 min; EI 70 eV, m/z 25-250. MS used for identity and for by-products.
- **Quantification sentence (verbatim):** "Binding percentage of flavours was determined from the
  difference between the peak areas of flavoured samples in the absence and presence of proteins
  such that: Binding% = (1 - Peak area with protein added / Peak area without protein added) x
  100%." This is headspace depletion: it counts reversible hydrophobic binding, covalent binding and
  any consumption by side reactions alike, and it assumes the protein does not change the
  air-liquid partition of the free flavour.
- **DSC (PPIs only).** 10 % w/v PPIs in 0.3 M NaCl (a different medium from the binding assay),
  flavour 100 / 250 / 500 ppm, 1 h rotary shaking; 10-15 uL in hermetic Tzero pans, 30 -> 120 C at
  10 C/min; delta-H and Td from the endotherm; duplicates.
- **Design.** Study 1: 2 x 3 x 2 factorial (class x carbon number x extraction) for canola and pea,
  gluten added for comparison, duplicated. Study 2 (DSC): 2 x 3 x 3 (class x carbon number x
  concentration). Tukey p < 0.05.

## 3. Tables re-typed

The paper has **no tables**. Its numerical results are Fig. 1a (percentage bound of aldehydes, five
proteins x C6/C7/C8), Fig. 1b (ketones), Fig. 5 (delta-H of PPIs with each flavour at 250 ppm) and
Fig. 6 (delta-H vs 100/250/500 ppm). All four are **FIGURE-ONLY**. What the text layer contains from
them is axis ticks and the Tukey superscript letters, not values; the aldehyde panel is drawn on a
0-80 % axis and the ketone panel on a 0-30 % axis, which is the only scale information available
and is not a measurement.

Numbers printed in the running text (all of them):

| item | as printed | where |
|---|---|---|
| Protein contents | CPIs 87.32 %, CPIa 75.35 %, PPIs 82.68 %, PPIa 82.82 %, gluten 76.01 % (N x 5.7) | §2.5 |
| Aldehyde vs ketone | "all proteins exerted higher binding capacities to aldehydes than ketones; binding capacities were 2-5 times higher (Fig. 1a vs. b)" | §3.1.2 |
| Aldehyde ordering | "CPIs revealed the highest binding capacity to aldehydes, followed by wheat gluten and PPIs" | §3.1.3 |
| Ketone ordering | "PPIs > wheat gluten > CPIs, although 2-hexanone was an exception and showed the lowest affinity to PPIs" | §3.1.3 |
| Chain length | "retention of both aldehydes and ketones was significantly enhanced with an increase in flavour carbon number" (all proteins) | §3.1.1 |
| Extraction method | "PPIa had a higher binding capacity to aldehydes compared with PPIs"; CPIa bound "much less" aldehyde than CPIs; ketones bound more to salt-extracted proteins except 2-octanone on CPIa | §3.1.4 |
| Pea amino acids (cited from Khattab, Arntfield & Nyachoti 2009) | Cys 0.35, Met 1.60, Lys 6.25 g/100 g protein (PPIs) | §3.1.3 |
| Canola 12S globulin (cited from Schwenke 1981) | Cys 1.07, Met 1.84, Lys 3.45 g/100 g | §3.1.3 |
| By-products (CPIs only) | 2-butyl-2-octenal (C12H22O) from hexanal, RT 11.77 min; 2-pentyl-2-nonenal (C14H26O) from heptanal, RT 12.72 min; "No volatile flavour by-products were detected in other protein-flavour mixtures"; "No such products were observed when octanal was the added aldehyde" | §3.1.5 |
| DSC | Td unchanged by flavours ("data not presented"); delta-H lower than control for every flavour, ketones higher delta-H than the matching aldehyde; delta-H falls linearly with 100 -> 500 ppm | §3.2 |
| Cited pea vicilin (Heng 2004) | 0.1 % vicilin: pentanal 75 % -> octanal 88 % retained; 2-octanone bound 16 % more than 2-pentanone | §3.1.1 |
| Cited soy (Gremli 1974) | 5 % soy protein binds aldehydes 2-4x more than ketones of the same carbon number | §3.1.2 |

## 4. Numbers the repository can use

| quantity | value | unit | conditions | source location | evidence class | `binding_constants.yml` fit |
|---|---|---|---|---|---|---|
| Hexanal / heptanal / octanal % bound, PPIs | FIGURE-ONLY | % | 10 g powder/L (8.27 g protein/L), 250 ppm v/v (~2.0 / 1.8 / 1.6 mM), 0.01 M K-phosphate pH 8, 30 C, 3 h, headspace at 40 C | Fig. 1a | figure_only | would be `percent_bound_at_conditions`, `protein_basis: g_isolate_powder`, `protein_purity_fraction: 0.8268`; **cannot be filled** |
| Same, PPIa, CPIs, CPIa, gluten | FIGURE-ONLY | % | as above (PPIa, CPIa, gluten as 1 % suspensions) | Fig. 1a | figure_only | as above |
| 2-hexanone / 2-heptanone / 2-octanone % bound, all five proteins | FIGURE-ONLY | % | as above | Fig. 1b | figure_only | as above |
| Aldehyde : ketone binding ratio, same carbon number | 2-5 | x | all proteins, conditions above | §3.1.2 text | measured (ratio, printed as a range) | no record type; a within-study ratio only |
| PPIs protein purity | 82.68 | % (N x 5.7) | freeze-dried powder | §2.5 | measured | already carried as `protein_purity_fraction: 0.8268` |
| PPIa protein purity | 82.82 | % (N x 5.7) | freeze-dried powder | §2.5 | measured | new; no PPIa record exists |
| Lys content of PPIs | 6.25 | g / 100 g protein | cited from Khattab et al. 2009, not measured here | §3.1.3 | level_only, secondary | not a binding record. Arithmetic for `protein_matrices.yml` (amine): 6.25 g / 146.19 g/mol = **0.428 mmol Lys per g protein** = 0.354 mmol per g PPIs powder (x 0.8268). Would need Khattab 2009 read before it is entered; it is total lysine, not free epsilon-amine. |
| Cys content of PPIs | 0.35 | g / 100 g protein | cited, as above | §3.1.3 | level_only, secondary | 0.35 / 121.16 = **0.0289 mmol half-cystine per g protein** (free SH + 2 x SS together; the split is not given here — Gao 2020 Table 2 gives the split for an alkaline isolate, see `gao2020_extraction.md`) |
| Met content of PPIs | 1.60 | g / 100 g protein | cited | §3.1.3 | level_only, secondary | none |
| Aldol by-products with canola isolate | 2-butyl-2-octenal (hexanal), 2-pentyl-2-nonenal (heptanal); none with PPIs, PPIa, gluten; none from octanal | identity only | 1 % CPIs, 250 ppm, pH 8, 30 C, 3 h | §3.1.5, Figs 2-4 | measured (qualitative) | none; it is mechanistic support for an amine (Schiff base) route being active at 30 C on a plant globulin |
| delta-H of PPIs with / without flavours | FIGURE-ONLY | J/g | 10 % PPIs, 0.3 M NaCl | Figs 5-6 | figure_only | none |

Directional claims the model can hold as hold-out shapes (no numbers): (i) binding rises with chain
length C6 < C7 < C8 for both classes on every protein; (ii) aldehydes bind 2-5x more than the
ketone of the same chain; (iii) for pea, the alkaline/acid isolate binds aldehydes MORE than the
salt isolate, the opposite of canola; (iv) aldehyde addition lowers the denaturation enthalpy of
pea isolate in proportion to concentration, ketones less so.

## 5. Flags

1. **Every percent-bound value is FIGURE-ONLY.** The paper prints no table and no numeric binding
   value in its text. The 2014 paper cannot supply a pea hexanal `percent_bound_at_conditions` row;
   this agrees with the FIGURE-ONLY warning already in `binding_constants.yml` for thesis Fig. 3.1a.
2. **Depletion, not affinity.** Binding % is a single-point headspace depletion at 40 C after 3 h at
   30 C, with no internal standard and no calibration. It cannot separate reversible from covalent
   binding, and the 40 C sampling step sits 10 C above the equilibration temperature.
3. **"ppm" is volumetric.** The 1000 ppm stock is 0.1 mL per 100 mL, so 250 ppm = 250 uL/L =
   203-205 mg/L depending on density (§2), not 250 mg/L. The molarity differs across the homologous
   series (2.03 mM for C6 down to 1.60 mM for C8), so the chain-length comparison is at equal volume,
   not equal moles.
4. **Two loading forms.** PPIs was a dispersed 2 % solution diluted to 1 %; PPIa and gluten were dry
   powder weighed into the vial (suspension). The PPIa-over-PPIs aldehyde result therefore also
   compares a dispersion with a suspension.
5. **Existing repo records.** `wang2015_ppi_*` records in `binding_constants.yml` carry
   `flavor_concentration_mg_per_L: 250.0`; the two journal papers read here prepared stocks
   volumetrically. Whether the thesis Chapter 5 stock was volumetric too was not checked in this
   pass (thesis not re-read); if it was, that field should read ~204 mg/L (2.0 mM for hexanal).
6. **Column film "20 um"** as printed for a 0.2 mm i.d. VF-5ms; 0.2 um is the plausible value.
   Irrelevant to the numbers, noted for fidelity.
7. **Lys / Cys / Met for pea are cited, not measured** (Khattab et al. 2009, LWT 42:1107). The amine
   density computed in §4 (0.428 mmol/g protein) is offered as arithmetic on a secondary number and
   must not enter `protein_matrices.yml` without reading the primary source. It is total lysine, and
   the engine's amine band for BLG already declares that not all lysines are free.
8. **The by-products argue for covalent chemistry on canola but were NOT seen on pea** at 30 C /
   3 h; the authors attribute pea's weaker aldehyde binding to its low Cys+Met (not to Lys, which is
   higher in pea than canola). For the engine, which routes hexanal to the amine pool, this is a
   caution: in this dataset the protein richer in amine bound less aldehyde.
9. **No pH, salt or temperature series** in this paper (all at pH 8, 0.01 M phosphate, 30 C); the
   pH/salt series is the FRI 2015 paper, the heat series is the Food Hydrocolloids 2015 paper.
10. **No time course** beyond the statement that 3 h sufficed; the preliminary data are not shown.
