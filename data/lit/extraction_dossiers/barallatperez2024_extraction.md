# Barallat-Perez 2024 — EXTRACTION (commercial lupin protein isolate 1 % w/v in MilliQ water pH 7.0, hexanal / nonanal / 2-nonanone at 5 mg/L; in vivo nose-space PTR-ToF-MS with simultaneous time-intensity sensory on ten subjects at 25 C, plus an in vitro static-headspace GC-MS binding assay at 30 C with and without 0.01 % pig gastric mucin)

### THE PAPER THE REPOSITORY ALREADY SEALED, READ AT LAST — and it turns out the seal was on the wrong document: `HOLDOUT_SEALED_BINDING` in `src/kinetic_core/parameters_matrix.py` names six lupin and mucin constants, and **this paper prints no binding constant and no binding percentage at all** (every binding number is inside Figure 4), so those sealed values must come from the group's earlier 2023 paper; what this paper does deliver is the corpus's only IN-MOUTH measurement — a matched pair of instrumental release and trained-panel perception on the same swallow — and lupin is not a protein `data/species/protein_matrices.yml` carries.

**Source on disk:** `data/articles/Barallat-Perez2024.pdf` (11 pp., J. Agric. Food Chem. 2024, 72,
8731-8741, CC-BY 4.0).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Barallat-Perez2024.txt`). **Tables 1, 2 and 3 are all printed SIDEWAYS on the
page and came through the text layer scrambled** — the values extracted as bare columns with their
row and column labels detached, and Table 1 is an image with no text layer at all. All three were
therefore **re-read from rendered page images** (pages 2, 5 and 7 at 170-200 dpi) and are re-typed
in full below from those renders; the values agree with the scrambled text-layer extraction
wherever the extraction produced anything at all. **Equations 1, 2 and 3 also lost their minus signs
in extraction** and were re-read from the page-4 render: eq 1 is `I = a t^(-b)`, not `I = a t^b`.
Figures 1 (study schematic), 2A-C, 3A-C (release and perception curves), 4 (the binding
percentages), S1, S2 and S3 are images. **Every binding percentage in this paper is figure-only.**
Supporting Information (Table S1, Table S2, Figures S1-S3) is **not on disk**.

Repo status before this dossier: Barallat-Perez appears in `data/lit/binding_constants.yml` under
two OTHER source ids — `barallat_perez_2023_jafc` (JAFC 2023, 71(50), 20274-20284) and
`barallat_perez_2025_npjsf` — and in `HOLDOUT_SEALED_BINDING` (`parameters_matrix.py`) as six
key-only entries. **This 2024 paper has no source id anywhere in the repository and no extraction
dossier.**

## 0. Identity

| field | value |
|---|---|
| Title | "Drivers of the In-Mouth Interaction between Lupin Protein Isolate and Selected Aroma Compounds: A Proton Transfer Reaction−Mass Spectrometry and Dynamic Time Intensity Analysis" |
| Authors | Cristina Barallat-Pérez (corresponding, cristina1.barallatperez@wur.nl, ORCID 0000-0003-0963-2242), Michele Pedrotti (Fondazione Edmund Mach, San Michele all'Adige), Teresa Oliviero, Sara Martins (also AFB International EU, Oss), Vincenzo Fogliano, Catrienus de Jong (Wageningen Food and Biobased Research) — Department of Agrotechnology and Food Science, Wageningen University & Research |
| Venue | J. Agric. Food Chem. 2024, 72 (16), 8731-8741. Received 24 November 2023, revised 15 March 2024, accepted 22 March 2024, published 5 April 2024. Part of the virtual special issue "13th Wartburg Symposium on Flavor Chemistry and Biology" |
| DOI | 10.1021/acs.jafc.3c08819 |
| Licence | **CC-BY 4.0** ("This article is licensed under CC-BY 4.0"), © 2024 The Authors, published by the American Chemical Society |
| The three ligands | **hexanal** (CAS 66-25-1, C6H12O), **nonanal** (CAS 124-19-6, C9H18O), **2-nonanone** (CAS 821-55-6, C9H18O) — a chain-length pair (hexanal vs nonanal) and a carbonyl-position pair (nonanal vs 2-nonanone) at constant C9 |
| Protein | **Lupin Protein Isolate 10600**, ProLupin GmbH, Grimmen, Germany. Manufacturer's specification: **91 % protein, 3 % lipid**; aqueous extraction and spray drying from seeds of sweet blue lupin, *Lupinus angustifolius* L.; taste from neutral (pH 7.0) to grassy, grainy and flour-like odour |
| Second protein | **pig gastric mucin** (Sigma-Aldrich), the salivary component, used in artificial saliva and in the in vitro assay |
| Naming | AUC_R / Imax_R / Tmax_R = the in vivo PTR-ToF-MS release parameters; AUC_S / Imax_S / Tmax_S = the same three on the trained sensory attribute "green"; `a` = fitted initial intensity and `b` = fitted decay rate from I = a t^(−b); "lingering" = persistence in seconds after the last swallow |
| Companion papers by the same group | ref 7 = **Barallat-Pérez, Janssen, Martins, Fogliano, Oliviero, JAFC 2023, 71(50), 20274-20284**, "Unraveling the Role of Flavor Structure and Physicochemical Properties in the Binding Phenomenon with Commercial Food Protein Isolates" — the in vitro parent study, and the one the repository's sealed lupin constants must actually come from |
| Companions on disk | `bi2022_extraction.md` (pea, the same headspace-binding family), `guo2020_extraction.md` (soy), `bornhorst2017_extraction.md` / `bornhorst2017b_extraction.md` (whey), `damodaran1981_extraction.md`, `andriot2000_extraction.md`, `leksrisompong2010_extraction.md`, `Meynier2002_extraction.md` |

## 1. Why it matters

**Everything measured here is NON-COVALENT retention, and the paper is explicit that the covalent
channel is a possibility it did not measure**: "Aldehydes can bind to proteins through reversible or
irreversible mechanisms, such as cysteine-aldehyde condensation reactions and Schiff base formation
under certain conditions (e.g., pH 6−10), forming strong amide linkages." The systems sit at pH 7.0,
inside that window, for 3 h at 30 C — so the aldehyde numbers are an **upper bound on reversible
retention with an uncontrolled covalent share**, exactly like Meynier's quarantined t-2-hexenal row.
The one clean control in the paper is **2-nonanone**, a ketone, which the authors say "predominantly
bind[s] through weaker hydrophobic interactions".

**Which layer.** These numbers belong to the matrix-retention side —
`src/kinetic_core/parameters_matrix.py` and the live measured registry
`data/lit/binding_constants.yml` (whose `percent_bound_at_conditions` record type is precisely the
shape of a Figure-4 binding percentage). **Nothing here goes to
`src/kinetic_core/matrix_sites.py`**: that module's `BINDING_CLASSES` are second-order covalent rate
constants in M^-1 s^-1 with activation-energy bands, and this paper measures no rate and no
temperature dependence of one.

**Is lupin a protein the matrix table carries? No.** `data/species/protein_matrices.yml` holds three
matrices — `blg` (beta-lactoglobulin from sequence counts), `soy_isolate` and `pea_isolate` (both
from measured thiol / disulfide / lysine densities). **There is no lupin entry and no mucin entry.**
`parameters_matrix.py`'s `MATRIX_LOADING` table likewise holds only `water`, `skim_milk`,
`caseinate_1pct`, `gelatin_3pct` and `soy_paste_hong`. So a lupin row from this paper would need a
matrix created first, and it would have **no site densities to pair with** — the covalent layer
would charge nothing for it and say so, which is the correct behaviour.

**The finding that most matters to the repository, and it is a negative one.**
`HOLDOUT_SEALED_BINDING` registers six keys — `kg_hexanal_lupin`, `kg_nonanal_lupin`,
`kg_2_nonanone_lupin`, `kg_hexanal_mucin`, `kg_nonanal_mucin`, `kg_2_nonanone_mucin` — each
annotated "Barallat-Perez 2024 lupin / pig gastric mucin — Module 6 STAR HOLD-OUT (D.6). Value not
carried in this file." Those three compounds are exactly this paper's three, and mucin is exactly
this paper's second protein, so the attribution is understandable. **But this paper prints no such
value.** Its nine binding percentages (M + aroma, LPI + aroma, LPI + M + aroma, for each of the
three compounds) exist only as bars in Figure 4 with no data labels, and the running text quotes
only a range ("increased 4−12 times following the addition of mucin"). A per-gram constant
`kg_hexanal_lupin` therefore cannot have been read from this document. Either it was read off
Figure 4's bars — which the house rules forbid and which `data/lit/binding_constants.yml`'s own
`pooling_caveat` on the 2023 record already refuses ("reading a bar chart is not a content-verified
number") — or it came from **ref 7, Barallat-Pérez et al. JAFC 2023**, which is the in vitro binding
paper and is already in `binding_constants.yml` as `barallat_perez_2023_jafc`. **The seal's citation
should be checked and, if it is the 2023 paper, corrected.**

**What this paper uniquely adds.** It is the only source in the five-paper cluster, and as far as
this dossier can see the only one in the corpus, that measures **release and perception on the same
swallow**. Every other binding study in `REVERSIBLE_BINDING` is an equilibrium headspace or dialysis
measurement on a sealed vial; this one puts the matrix in a mouth. Three consequences bear directly
on the layer's declared policy:

1. **The instrument and the panel disagree, and the paper says so.** Protein addition cuts nonanal's
   instrumental Imax_R by **72.41 %** while cutting the panel's Imax_S by only **15.23 %**. That is
   a factor of **4.8 (mine)** between the physical suppression and the perceived suppression on the
   same samples. `parameters_matrix.py` reports ratios and rankings rather than absolute ppb
   precisely because the corpus's ratios span 2000x; this paper adds the further warning that
   **even a correct headspace suppression does not transfer to a perceived one**.
2. **The unflavoured lupin blank already scores "green" at 51 ± 5 on a 100 mm scale** (Table 2), a
   quarter below flavoured hexanal's 67 ± 4. The sensory signal is riding on a large matrix baseline,
   and hexanal and nonanal "were not detected in unflavored samples in vivo" — so the panel is
   scoring lupin's own off-notes into the same attribute. Any threshold or OAV work fed from this
   paper inherits that.
3. **Saliva is not neutral.** Adding 0.01 % w/v pig gastric mucin to the in vitro system raises
   binding 4-12x, and the authors observe that the LPI + mucin binding is **greater than the sum of
   the two separately**. The matrix layer has no salivary term at all, and this is a measured reason
   why a vial-derived retention will under-predict an in-mouth one.

What this paper does NOT give the repository: any binding constant in any unit; any printed binding
percentage; any protein molar mass; any temperature series; any pH other than 7.0; any covalent
adduct measurement; any odour threshold; any concentration in ppb of an aroma compound in the liquid
phase; any lupin site density.

## 2. Methods as they matter to a model

- **The protein and its exact description.** Lupin Protein Isolate 10600 (ProLupin GmbH, Grimmen,
  Germany), a **commercial** isolate — not laboratory-made, unlike Bi 2022's pea. Manufacturer's
  specification: **91 % protein and 3 % lipid**, obtained by **aqueous extraction and spray drying**
  from sweet blue lupin (*Lupinus angustifolius* L.). Batches stored at **10-15 C**, dry, away from
  light and air. **This is the only isolate in the cluster whose protein content is stated**, which
  matters because `binding_constants.yml`'s units discipline distinguishes `g_protein` from
  `g_isolate_powder` and calls mixing them an up-to-1.3x error: here the conversion factor is known
  and is **0.91**.
- **Stock and loading.** Stock at **2 % w/v in MilliQ water, pH 7.0** — *water, not buffer*;
  vortexed 10-20 s at 3200 rpm, held **20 min at 30 C** in a water bath, vortexed again. The model
  systems are **0 or 1 % w/v LPI**, i.e. **10 g of isolate powder per litre = 9.1 g protein per
  litre (mine)**.
- **Aroma.** Hexanal, nonanal and 2-nonanone (Sigma-Aldrich, Zwijndrecht), purity >= 95 %. Stocks at
  **10 mg/L in MilliQ water pH 7.0** in 100 mL amber glass, held **1 h at 30 C**. **Final
  concentration in every system 5 mg/L**, each compound added separately (never as a mixture),
  chosen as below the FEMA GRAS 25th-edition maximum use level. **No co-solvent is mentioned**,
  which distinguishes this preparation from Bi 2022's ~1.25 % methanol.
- **The seven systems.** Three aroma-only (no protein), three aroma + 1 % LPI, one LPI-only blank.
  **10 mL** of aqueous model system, incubated in a **shaking water bath at 125 rpm for 3 h** before
  nose-space analysis; 3 h is stated to be "adequate timing for achieving equilibrium" on the
  authority of Wang & Arntfield.
- **Artificial saliva.** Adapted from van Ruth. Per 1000 mL: NaHCO3 5.208 g, K2HPO4·3H2O 1.369 g,
  NaCl 0.877 g, KCl 0.477 g, CaCl2·2H2O 0.441 g, **pig gastric mucin 2.160 g**, NaN3 0.5 g. The text
  also says "Artificial saliva was made at 0.01 wv%" and the GC-MS assay used "0.01 wv% mucin"
  (= 0.1 g/L). **2.160 g/L is 0.216 % w/v, twenty-one times the stated 0.01 % w/v** — the two
  statements cannot both describe the mucin concentration in the assay (Flags 4).
- **Method A — in vivo nose-space PTR-ToF-MS (what it measures).** High-sensitivity **PTR-QiToF-MS**
  (Ionicon Analytik, Innsbruck). Drift tube **100 C, 900 V, 460 Pa, E/N = 133 Td**. Sampling through
  a PEEK capillary (1/16" OD, 0.01" ID) heated to **100 C at 40 mL/min**. Mass resolution
  m/Δm >= 4800. Two Teflon nose pieces (6.8 mm diameter, 6.4 cm long) in the nostrils, connected to
  a heated (100 C) N.A.SE device. 20 s background, then 1 min of regular breathing for a baseline.
  **This measures the concentration of the compound in the exhaled nasal air during and after
  swallowing, in ppbV** — a non-equilibrium, saliva-diluted, body-temperature quantity, not a
  partition coefficient.
- **Mass channels.** Hexanal quantified at **m/z 101.103** with primary fragment **m/z 83.055**;
  nonanal and 2-nonanone both at **m/z 143.158** with primary fragment **m/z 125.142**. "Absolute
  quantification was derived by summing the obtained values corresponding to the molecular ion
  fragments." **Nonanal and 2-nonanone share a channel and are therefore never measured in the same
  sample** — each was dosed separately, which is why the design has one compound per system.
- **Method B — dynamic sensory time intensity (what it measures).** EyeQuestion v5. **Ten European
  female subjects, 26 ± 2 years**, non-smoking, no swallowing disorder, no lupin allergy, no dental
  braces; **saliva flow rate 0.145 ± 0.1 g/min** and **mouth volume 75 ± 8.5 g water** measured.
  Three training sessions; the single trained attribute is **"green"**, defined as "reminiscent of
  grass and vegetables, with a slight pungency, accompanied by hints of fruitiness and freshness".
  Samples served at **25 ± 5 C** in a 20 mL clear GC-MS glass vial (75.5 x 17.5 mm), sipped through
  a straw, **held in the mouth 10 s, swallowed, a second swallow 10 s later, occasionally a third**.
  Rated on a **100 mm unstructured line scale** anchored "very weak" to "very strong". Maximum six
  samples per session; palate cleansed with water and unsalted crackers. **Evaluated in triplicate,
  n = 10 subjects.** Study exempted from ethical approval by the Wageningen medical ethics
  committee; Declaration of Helsinki; written consent under EU 679/2016; participants paid.
- **Panel recruitment.** A prior focus group of **n = 40** regular plant-based-beverage consumers at
  Wageningen chose between soy, lupin and pea isolate on overall taste and odour; **LPI won with
  52.5 %** (Figure S1). So lupin's selection is itself a measured preference, not an arbitrary
  choice.
- **Method C — in vitro static headspace GC-MS binding (the only equilibrium measurement here).**
  Agilent 7890A GC coupled to an Agilent 5975C MSD with triple-axis detector. Three sample
  combinations — **LPI + aroma**, **mucin + aroma**, **LPI + aroma + mucin** — at **1 % w/v LPI,
  5 mg/L aroma, 0.01 % w/v mucin**, plus buffered LPI and mucin references with no aroma. Sealed
  vials, **water-bath shaker at 30 C and 125 rpm for 3 h**, **triplicate**. Binding is computed by
  headspace depletion against the protein-free control:

  > eq 2: `binding(%) = (1 − (HS1 − HS2)/HS3) × 100`
  > eq 3: `binding(%) = 1 − ((HS1 − HS2 − HS4 − HS5)/HS3) × 100`

  where HS1 is the flavoured protein solution's headspace abundance, HS2 the headspace without
  aroma, HS3 the headspace without protein, HS4 the mucin-solution-plus-buffer headspace without
  aroma, and HS5 "the headspace abundance of the protein-based mucin solution". **Both equations are
  re-typed here from the page image; the minus signs did not survive text extraction.** Eq 3's outer
  bracket is missing in the printed article, and its subtraction of both HS4 and HS5 is difficult to
  read as a blank correction (Flags 5).
- **Method D — lingering and decay.** Curves fitted to **I = a t^(−b)** (eq 1), `a` the initial
  intensity and `b` the decay rate, applied separately to the PTR-ToF-MS release curve and to the
  sensory curve. Lingering is averaged **per second, over all subjects and replicates, from after
  the third and last swallow to the end of the test**, and is reported in seconds.
- **Statistics.** Two-way ANOVA per sample combination on AUC, Imax and Tmax, with Tukey post hoc at
  **p < 0.05**; GraphPad Prism 9.3.1471 and RStudio 4.2.1. Curves smoothed with `geom_smooth` in
  ggplot2 — **the plotted curves are smoothed, the tabulated parameters are not**.

## 3. Tables re-typed

All three tables are printed sideways and were read from rendered page images (see the source note).

### Table 1 (p. 8732). "Physicochemical and Structural Features of the Selected Aroma Compounds"

Footnote a: "(1−5) Properties obtained from ref 31" — **ref 31 is PubChem**, accessed 2024-02-09.
So no property in this table was measured by the authors.

| Selected aroma compounds | CAS | Molecular weight (g/mol) | Vapor pressure (mmHg) | LogP | Water solubility (mg/L) |
|---|---|---:|---:|---:|---:|
| Hexanal (C6H12O) | 66-25-1 | 100.16 | 11.26 | 1.8 | 5640 |
| Nonanal (C9H18O) | 124-19-6 | 142.24 | 0.37 | 3.3 | 96 |
| 2-nonanone (C9H18O) | 821-55-6 | 142.24 | 0.62 | 3.1 | 371 |

(The "Chemical structure" column is a drawn skeleton and carries no number.)

### Table 2 (p. 8735). "Summary of Parameters (Mean ± SE) Describing the In Vivo Hexanal, In Vivo Nonanal, In Vivo 2-Nonanone, and Dynamic 'Green' Perceived Intensity for Flavored Lupin Protein-Based Aqueous Model Systems"

Footnote a: "Letters denote significant differences (p < 0.05). Treatments with the same letter are
not significantly different."

The LPI blank was measured on **both** mass channels; the flavoured samples on the channel of their
own compound. Release parameters (AUC_R, Imax_R, Tmax_R) are in ppbV-based units and seconds; the
paper prints no unit on AUC or Imax. Sensory parameters (AUC_S, Imax_S, Tmax_S) are on the 100 mm
scale and in seconds.

| parameter | LPI (m/z 101.103) | LPI (m/z 143.158) | hexanal (m/z 101.103) | LPI + hexanal (m/z 101.103) | nonanal (m/z 143.158) | LPI + nonanal (m/z 143.158) | 2-nonanone (m/z 143.158) | LPI + 2-nonanone (m/z 143.158) |
|---|---|---|---|---|---|---|---|---|
| AUC_R | 251 ± 11 (d) | 243 ± 10 (d) | 727 ± 69 (b) | 581 ± 41 (bc) | 401 ± 36 (cd) | 271 ± 11 (d) | 1515 ± 116 (a) | 1261 ± 116 (a) |
| Imax_R | 6 ± 2 (d) | 5 ± 1 (d) | 207 ± 32 (b) | 142 ± 18 (bc) | 58 ± 13 (cd) | 16 ± 2 (d) | 209 ± 28 (a) | 132 ± 15 (a) |
| Tmax_R | 12 ± 2 (b) | 26 ± 3 (a) | 6 ± 1 (b) | 5 ± 4 (b) | 10 ± 2 (b) | 7 ± 1 (b) | 9 ± 1 (b) | 8 ± 1 (b) |
| AUC_S | 1800 ± 268 (ab) | — | 2099 ± 222 (ab) | 1555 ± 211 (b) | 2735 ± 248 (a) | 2127 ± 287 (ab) | 2245 ± 209 (ab) | 2155 ± 295 (ab) |
| Imax_S | 51 ± 5 (b) | — | 67 ± 4 (ab) | 54 ± 5 (b) | 71 ± 4 (a) | 60 ± 5 (ab) | 70 ± 4 (ab) | 62 ± 5 (ab) |
| Tmax_S | 6 ± 2 (a) | — | 7 ± 1 (a) | 9 ± 2 (a) | 10 ± 3 (a) | 10 ± 3 (a) | 6 ± 1 (a) | 7 ± 2 (a) |

The sensory rows are single-valued for the LPI blank (there is one sensory measurement of it, not
one per mass channel); the printed table places the blank's sensory letters in the m/z 143.158 letter
column, and they are transcribed above against the m/z 101.103 sub-column where the values sit.

### Table 3 (p. 8737). "Initial Intensity (a), Decay Rate (b), Lingering Duration for All Samples and Subjects, In Vivo Aroma Release (PTR-ToF-MS_R), and Sensory Perception (Sensory_S) (Eq 1)"

Footnote a: "Data are presented as the average of the three replicates with the standard error.
Letters denote significant differences (p < 0.05). Treatments with the same letter are not
significantly different." **The `a` and `b` rows carry no standard error and no significance
letter**; only the lingering row is lettered.

| row | hexanal | LPI + hexanal | nonanal | LPI + nonanal | 2-nonanone | LPI + 2-nonanone |
|---|---:|---:|---:|---:|---:|---:|
| PTR-ToF-MS_R, a | 4.252 | 3.925 | 3.990 | 3.654 | 11.138 | 9.509 |
| PTR-ToF-MS_R, b | 0.002 | 0.001 | 0.002 | 0.002 | 0.012 | 0.008 |
| Sensory_S, a | 11.426 | 8.031 | 18.985 | 14.238 | 7.830 | 7.830 |
| Sensory_S, b | 0.016 | 0.016 | 0.015 | 0.010 | 0.018 | 0.018 |
| lingering (s) | 60 ± 10 (ab) | 71 ± 10 (c) | 88 ± 6 (a) | 75 ± 8 (ab) | 74 ± 9 (ab) | 75 ± 9 (ab) |

Two oddities in the printed table, transcribed as they stand: **2-nonanone and LPI + 2-nonanone
carry identical Sensory_S parameters (a = 7.830, b = 0.018)**, and hexanal / LPI + hexanal carry
identical Sensory_S b (0.016), as do 2-nonanone's pair. The lingering significance letters are also
inconsistent with the values as printed (LPI + hexanal at 71 s is lettered "c" while hexanal at 60 s
is "ab"), which cannot be read as a simple ordering (Flags 6).

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| chain length, hexanal -> nonanal: AUC_R decrease | **44.89 %** | Results, "Chain Length" |
| chain length, hexanal -> nonanal: Imax_R decrease | **71.92 %** | same, and the abstract |
| Tmax_R over chain length | no significant differences | same |
| protein addition, AUC_R decrease: LPI + hexanal / LPI + nonanal | **20.06 % / 32.37 %** | same |
| protein addition, Imax_R decrease: LPI + hexanal / LPI + nonanal | **30.91 % / 72.41 %** | same, and the abstract |
| protein-free chain length, AUC_S / Imax_S increase | **30.31 % / 6.24 %** | same |
| protein addition, AUC_S decrease (hexanal / nonanal) | **25.91 % / 22.25 %** | same |
| protein addition, Imax_S decrease (hexanal / nonanal) | **25.92 % / 15.23 %** | same (but see the arithmetic below — the hexanal figure does not reproduce) |
| carbonyl position, 2-nonanone -> nonanal: AUC_R reduction | **73.52 %** | Results, "Reactivity and Position of the Carbonyl Group" |
| carbonyl position, 2-nonanone -> nonanal: Imax_R reduction | **72.25 %** | same, and the abstract |
| protein addition to 2-nonanone: AUC_R / Imax_R decrease | **16.83 % / 36.84 %** | same, and the abstract |
| protein-free carbonyl position, AUC_S / Imax_S increase | **21.8 % / 1.87 %** | same |
| protein addition to 2-nonanone: AUC_S / Imax_S decrease | **4 % / 11.05 %** | same |
| hexanal volatility vs nonanal | "thirty-fold higher" | Results (from Table 1: 11.26/0.37 = 30.4, mine) |
| hexanal lingering vs nonanal | hexanal **46.94 % less persistent** than nonanal | Results, "Effect of Aroma Physicochemical Properties" |
| nonanal lingering vs 2-nonanone | nonanal exceeds 2-nonanone by **15.93 %** | same |
| effect of mucin on the in vitro binding response | **increased 4−12 times** | Results, "Effect of Mucin"; and the abstract |
| mucin + protein synergy | "the resulting binding effect does not simply sum up equally and proportionally ... higher binding than what would be expected solely on the basis of their individual contributions" | same |
| mucin loading rationale | "a minimal amount of mucin (0.01 wv%)" chosen because oral mucin varies with age, oral health and genetics | same |

**BINDING PERCENTAGES ARE FIGURE-ONLY.** Figure 4 shows nine bars — M + hexanal, LPI + hexanal,
LPI + M + hexanal; the same three for nonanal; the same three for 2-nonanone — with error bars and
significance letters (c, b, a / bc, a, a / c, b, a) and **no data labels**. Per house rule they are
not typed as numbers. The release and perception curves (Figures 2A-C, 3A-C) are likewise
figure-only, and they are smoothed.

### Arithmetic on the printed numbers (all mine)

**1. The printed percentage changes reproduce from Table 2, with one exception.** Using the printed
means:

| claim | printed | recomputed from Table 2 (mine) | verdict |
|---|---:|---:|---|
| AUC_R, hexanal -> nonanal | 44.89 % | (727−401)/727 = 44.84 % | agrees |
| Imax_R, hexanal -> nonanal | 71.92 % | (207−58)/207 = 71.98 % | agrees |
| AUC_R, protein on hexanal | 20.06 % | (727−581)/727 = 20.08 % | agrees |
| AUC_R, protein on nonanal | 32.37 % | (401−271)/401 = 32.42 % | agrees |
| Imax_R, protein on hexanal | 30.91 % | (207−142)/207 = **31.40 %** | 0.5 pt off (rounding of the means) |
| Imax_R, protein on nonanal | 72.41 % | (58−16)/58 = 72.41 % | exact |
| AUC_R, 2-nonanone -> nonanal | 73.52 % | (1515−401)/1515 = 73.53 % | agrees |
| Imax_R, 2-nonanone -> nonanal | 72.25 % | (209−58)/209 = 72.25 % | exact |
| AUC_R, protein on 2-nonanone | 16.83 % | (1515−1261)/1515 = 16.77 % | agrees |
| Imax_R, protein on 2-nonanone | 36.84 % | (209−132)/209 = 36.84 % | exact |
| AUC_S, protein on hexanal | 25.91 % | (2099−1555)/2099 = 25.92 % | agrees |
| AUC_S, protein on nonanal | 22.25 % | (2735−2127)/2735 = 22.23 % | agrees |
| Imax_S, protein on nonanal | 15.23 % | (71−60)/71 = 15.49 % | agrees |
| **Imax_S, protein on hexanal** | **25.92 %** | **(67−54)/67 = 19.40 %** | **DOES NOT REPRODUCE (Flags 3)** |

The failing cell's printed value, 25.92 %, is numerically the AUC_S figure 25.91 % to two decimals,
which is the signature of a copy error in the manuscript rather than a different computation.

**2. The lingering percentages use two different denominators.** "hexanal was 46.94 % less
persistent than nonanal": (88−60)/60 = **46.67 %**, i.e. relative to *hexanal*. "nonanal ...
surpassing 2-nonanone by 15.93 %": (88−74)/88 = **15.91 %**, i.e. relative to *nonanal*. Both
reproduce, but on opposite bases, so the two figures are not comparable to each other.

**3. Instrument against panel, the same samples.** Ratio of the instrumental Imax_R suppression to
the sensory Imax_S suppression on protein addition:

| compound | Imax_R suppression | Imax_S suppression (recomputed) | ratio (mine) |
|---|---:|---:|---:|
| hexanal | 31.4 % | 19.4 % | **1.6x** |
| nonanal | 72.4 % | 15.5 % | **4.7x** |
| 2-nonanone | 36.8 % | 11.4 % (= (70−62)/70) | **3.2x** |

**The panel is between 1.6 and 4.7 times less sensitive to the protein's suppression than the
instrument is.** This is the number a matrix layer that reports "what a sensory panel would notice"
must confront, and it is the clearest measurement of it in the cluster.

**4. Converting Table 2 to the registry's per-gram form — and why it should not be shipped.** The
registry's K_g = (K_water/K_matrix − 1)/protein_g_per_L. If one reads the with-protein / no-protein
release ratio as a headspace-suppression ratio at 10 g/L of isolate powder (or 9.1 g/L of protein):

| compound | Imax_R ratio | implied K_g, L per g powder | implied K_g, L per g protein | on AUC_R instead |
|---|---:|---:|---:|---:|
| hexanal | 142/207 = 0.686 | 4.6e-2 | 5.0e-2 | 2.5e-2 / 2.8e-2 |
| nonanal | 16/58 = 0.276 | 2.6e-1 | 2.9e-1 | 4.8e-2 / 5.3e-2 |
| 2-nonanone | 132/209 = 0.632 | 5.8e-2 | 6.4e-2 | 2.0e-2 / 2.2e-2 |

**These are NOT partition-derived constants and must not enter `REVERSIBLE_BINDING`.** A nose-space
intensity is measured after the sample has been diluted by saliva of unrecorded volume, warmed from
25 C to body temperature, aerated by breathing and partly swallowed; the AUC and Imax routes here
already disagree by up to **5.4x** on the same compound (nonanal), which is by itself proof that no
single equilibrium constant is being measured. The arithmetic is recorded only so that a later
reader can see it was done and rejected. Note also that the same compound gives **2.6e-1 L/g on
Imax_R and 4.8e-2 L/g on AUC_R** — the choice of release parameter moves the answer more than the
choice of protein does.

**5. The vapour-pressure ratio the text calls "thirty-fold".** 11.26/0.37 = **30.4** (hexanal over
nonanal) and 11.26/0.62 = **18.2** (hexanal over 2-nonanone). Confirmed against Table 1.

**6. The decay rates say the protein slows release, not perception.** From Table 3, the PTR decay
rate `b` falls on protein addition for hexanal (0.002 -> 0.001) and 2-nonanone (0.012 -> 0.008) and
is unchanged for nonanal (0.002 -> 0.002); the sensory `b` falls only for nonanal (0.015 -> 0.010)
and is unchanged for the other two. **The two channels do not even agree on which compound the
protein slows.** The initial intensity `a` falls on protein addition in every one of the six
comparisons except the identical 2-nonanone sensory pair.

**7. Lingering changes sign with the compound.** Protein addition raises hexanal's lingering
60 -> 71 s (**+18 %, mine**) and lowers nonanal's 88 -> 75 s (**−15 %, mine**), leaving
2-nonanone's essentially unchanged (74 -> 75 s). A matrix term that can only suppress is refuted by
the hexanal cell, in the same way `kg_delta_decalactone_caseinate` and `kg_furaneol_caseinate`
refute it in `REVERSIBLE_BINDING` — except that here the standard errors (± 10 and ± 8 s) make the
change marginal.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** `hexanal` is keyed. `nonanal` is keyed.
**2-nonanone is NOT in `compounds.yml`** — but it IS in `parameters_matrix.py`'s
`COMPOUND_STRUCTURE` as `2_nonanone` (methyl ketone, C9) and it already carries two shipped
constants, `kg_2_nonanone_soy` (3.72e-2 L/g, dialysis) and `kg_2_nonanone_blg` (6.63e-2 L/g,
headspace depletion), so the compound is live in the matrix layer while absent from the product
registry. **Lupin has no matrix key anywhere**; mucin has none either.

Every row below shares: commercial lupin protein isolate (ProLupin 10600, **91 % protein**) at
**1 % w/v = 10 g powder/L = 9.1 g protein/L**, aroma at **5 mg/L**, **MilliQ water, pH 7.0, no
buffer, no co-solvent**, 10 mL system, 3 h equilibration at 30 C with shaking at 125 rpm.
In vivo rows: served at **25 ± 5 C**, ten female subjects, triplicate. In vitro rows: sealed vial at
**30 C**, triplicate.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| in vivo release Imax_R, hexanal / LPI + hexanal | **207 ± 32 / 142 ± 18** | PTR-ToF-MS intensity (unit not printed; from ppbV) | 25 C in-mouth, m/z 101.103 | Table 2, p. 8735 | **measured_rate**-adjacent instrumental level; treat as `level_only` — it is not a rate and not a constant |
| in vivo release Imax_R, nonanal / LPI + nonanal | **58 ± 13 / 16 ± 2** | as above | m/z 143.158 | Table 2 | level_only |
| in vivo release Imax_R, 2-nonanone / LPI + 2-nonanone | **209 ± 28 / 132 ± 15** | as above | m/z 143.158 | Table 2 | level_only |
| in vivo release AUC_R, the six flavoured systems | **727 ± 69 / 581 ± 41 / 401 ± 36 / 271 ± 11 / 1515 ± 116 / 1261 ± 116** | as above | as above | Table 2 | level_only |
| in vivo release AUC_R / Imax_R / Tmax_R, LPI blank | **251 ± 11 / 6 ± 2 / 12 ± 2** (m/z 101.103) and **243 ± 10 / 5 ± 1 / 26 ± 3** (m/z 143.158) | as above, and s for Tmax | the unflavoured control | Table 2 | level_only — **the blank is not zero on either channel** |
| **Imax_R suppression by 1 % LPI: hexanal / nonanal / 2-nonanone** | **30.91 / 72.41 / 36.84** | % | 25 C in-mouth, 10 g powder/L | Results text (all three), and the abstract | **retention_percent** — but an IN-MOUTH, non-equilibrium suppression, not a headspace retention |
| **AUC_R suppression by 1 % LPI: hexanal / nonanal / 2-nonanone** | **20.06 / 32.37 / 16.83** | % | as above | Results text | **retention_percent**, same caveat |
| sensory Imax_S, the six flavoured systems and the blank | **67 ± 4 / 54 ± 5 / 71 ± 4 / 60 ± 5 / 70 ± 4 / 62 ± 5**, blank **51 ± 5** | mm on a 100 mm unstructured line scale, attribute "green" | 25 C, n = 10 in triplicate | Table 2 | **sensory** |
| sensory AUC_S, the six flavoured systems and the blank | **2099 ± 222 / 1555 ± 211 / 2735 ± 248 / 2127 ± 287 / 2245 ± 209 / 2155 ± 295**, blank **1800 ± 268** | mm·s | as above | Table 2 | **sensory** |
| sensory Tmax_S, all seven | **7 ± 1 / 9 ± 2 / 10 ± 3 / 10 ± 3 / 6 ± 1 / 7 ± 2**, blank **6 ± 2** | s | as above | Table 2 | **sensory** — no significant differences (all "a") |
| Imax_S suppression by 1 % LPI: hexanal / nonanal | **25.92 (does not reproduce; 19.40 mine) / 15.23** | % | as above | Results text | **sensory** — carry the recomputed 19.40 % for hexanal, or refuse the cell (Flags 3) |
| Imax_S suppression by 1 % LPI: 2-nonanone | **11.05** | % | as above | Results text | sensory |
| decay rate b, PTR release, six systems | **0.002 / 0.001 / 0.002 / 0.002 / 0.012 / 0.008** | dimensionless exponent in I = a t^(−b), t in s | 25 C in-mouth | Table 3, p. 8737 | **measured_rate** (a fitted decay exponent, not a first-order rate constant) |
| decay rate b, sensory, six systems | **0.016 / 0.016 / 0.015 / 0.010 / 0.018 / 0.018** | as above | as above | Table 3 | sensory |
| initial intensity a, PTR / sensory, six systems | PTR **4.252 / 3.925 / 3.990 / 3.654 / 11.138 / 9.509**; sensory **11.426 / 8.031 / 18.985 / 14.238 / 7.830 / 7.830** | fitted intercept of I = a t^(−b) | as above | Table 3 | level_only |
| aroma lingering, six systems | **60 ± 10 / 71 ± 10 / 88 ± 6 / 75 ± 8 / 74 ± 9 / 75 ± 9** | s after the last swallow | 25 C, n = 10 in triplicate | Table 3 | **measured_rate**-adjacent persistence; treat as level_only |
| **binding percentage, LPI + aroma, mucin + aroma, LPI + mucin + aroma (nine values)** | — | % | 30 C, 3 h, 1 % LPI, 5 mg/L aroma, 0.01 % mucin, static headspace GC-MS | **Figure 4, p. 8738** | **figure_only** — the nine bars carry no data labels |
| mucin's multiplicative effect on the in vitro binding | **4 to 12** | fold | as above | Results text and abstract | **within_study_ratio** — the only quantitative statement about Figure 4 that is printed |
| vapour pressure, hexanal / nonanal / 2-nonanone | **11.26 / 0.37 / 0.62** | mmHg | ambient (PubChem) | Table 1, p. 8732 | **level_only** — **NOT measured here**; from PubChem, ref 31 |
| LogP, hexanal / nonanal / 2-nonanone | **1.8 / 3.3 / 3.1** | — | PubChem | Table 1 | level_only, borrowed — **and the layer refuses any log P term (k4b guard #4)** |
| water solubility, hexanal / nonanal / 2-nonanone | **5640 / 96 / 371** | mg/L | PubChem | Table 1 | level_only, borrowed |
| lupin isolate protein content | **91** | % (with 3 % lipid) | manufacturer's specification | Materials, p. 8732 | level_only — **the only stated protein basis in the cluster**; gives `g_protein` = 0.91 x `g_isolate_powder` |
| instrument-vs-panel suppression ratio | **1.6x / 4.7x / 3.2x** for hexanal / nonanal / 2-nonanone | — | same samples, same swallow | (mine) from Table 2 | **within_study_ratio** — the headline transferable finding |
| chain-length contrast on Imax_R (C6 -> C9 aldehyde) | **71.92** | % reduction | protein-free, in-mouth | Results text | within_study_ratio |
| carbonyl-position contrast on Imax_R (2-nonanone -> nonanal) | **72.25** | % reduction | protein-free, in-mouth | Results text | within_study_ratio |
| per-gram binding constant on any basis | — | L/g | — | (mine, §3 item 4) | **derived_assumption, DO NOT SHIP** — the in vivo route gives 2.0e-2 to 2.9e-1 L/g depending on which release parameter is used, a 14x spread on the same data |

### Can these be put on the same basis as the shipped binding constants? Mostly not.

**(a) There is no equilibrium constant here to transport.** The only equilibrium measurement in the
paper is the in vitro GC-MS binding of Figure 4, and it is figure-only. The 4-12x mucin factor is a
ratio of two unreadable numbers and cannot be turned into a constant.

**(b) The in vivo suppressions are the wrong kind of object.** `REVERSIBLE_BINDING` multiplies a
per-gram constant by a protein loading to shift a headspace partition. A nose-space Imax is not a
headspace partition: saliva of unknown volume has diluted the matrix, the sample has warmed from
25 C toward 37 C, and breathing has stripped the mouth's air. The AUC and Imax routes disagree by up
to 5.4x on the same sample. **Carry these as observations of an in-mouth effect, in percent, with
their conditions — the shape `binding_constants.yml` calls `percent_bound_at_conditions` — and do
not invert them to a K_eff**, because the loading that the compound actually saw is not the 10 g/L
that was poured.

**(c) The one number that does transport is the protein basis.** 91 % protein is a manufacturer
specification for a named commercial isolate, and `binding_constants.yml`'s header explicitly wants
it: "a commercial isolate is 77-93 wt% protein, so mixing them is a up-to-1.3x error in c_p". ProLupin
10600 at 91 % sits at the top of that band.

**(d) Temperature.** 25 ± 5 C in the mouth (warming toward body temperature), 30 C in the vial. Both
are inside the existing spread of the shipped rows (25-40 C). Nothing here says anything about a
cooking temperature.

**(e) Nothing goes to `matrix_sites.py`.** No rate constant, no activation energy, no site density,
no thiol or amine assay on the lupin isolate.

## 5. Flags

1. **All three tables are printed sideways and the text layer scrambles them.** Table 1 has no text
   layer at all (it is an image). Tables 2 and 3 extract as detached columns of numbers with their
   row and column headers lost, and in Table 2 the letters and the values interleave. **Every value
   in section 3 was re-read from a rendered page image**; a future reader who trusts
   `pdftotext -layout` on this file alone will mis-assign cells. The equations lose their minus
   signs the same way — **eq 1 extracts as `I = at b` when the printed form is `I = a t^(−b)`**,
   which inverts the meaning of the decay parameter.
2. **The repository's sealed lupin and mucin constants are attributed to this paper and cannot come
   from it.** `HOLDOUT_SEALED_BINDING` in `parameters_matrix.py` names `kg_hexanal_lupin`,
   `kg_nonanal_lupin`, `kg_2_nonanone_lupin` and the three mucin twins as "Barallat-Perez 2024".
   **This paper contains no binding constant and no printed binding percentage.** The likely true
   source is ref 7, Barallat-Pérez et al., JAFC 2023, 71(50), 20274-20284, which is already carried
   in `data/lit/binding_constants.yml` as `barallat_perez_2023_jafc` — and whose own record there
   carries the caveat that its per-protein values "exist only in Figure 2 and are NOT transcribed
   here (reading a bar chart is not a content-verified number)". **Check the seal's citation, and
   check whether the sealed values were ever printed anywhere.**
3. **One printed percentage does not reproduce from the paper's own table.** "Imax_S decreased by
   25.92 %" for LPI + hexanal, against **19.40 %** recomputed from 67 ± 4 and 54 ± 5. The printed
   figure duplicates the AUC_S figure (25.91 %) to two decimals, which reads as a manuscript copy
   error. Every other percentage in the paper reproduces to within 0.5 points.
4. **The mucin concentration is stated two incompatible ways.** The artificial-saliva recipe puts
   **2.160 g of pig gastric mucin per 1000 mL (0.216 % w/v)**; the text says "Artificial saliva was
   made at 0.01 wv%" and the GC-MS assay used "0.01 wv% mucin" (0.1 g/L). These differ by **21.6x**.
   If the saliva is diluted to 0.01 % of *itself* the arithmetic does not work either. **Any K_eff
   inverted from a mucin binding percentage would be wrong by up to 21.6x**, and the binding
   percentages are figure-only in any case.
5. **Equation 3 is malformed as printed.** `binding(%) = 1 − ((HS1 − HS2 − HS4 − HS5)/HS3) × 100`
   has no outer bracket (so the ×100 does not multiply the whole expression as written), and
   subtracting **both** HS4 (mucin blank) and HS5 (protein + mucin blank) double-counts the mucin
   background. The intent is presumably a two-blank correction; as printed it is not
   dimensionally sensible against eq 2. Ask the authors, and do not re-derive any Figure 4 value
   from it.
6. **Table 3's significance letters do not order with its values, and two rows are duplicates.**
   LPI + hexanal at 71 ± 10 s lingering is lettered "c" while hexanal at 60 ± 10 s is "ab" and
   everything else is "ab" or "a" — so the letter that marks the *lowest* group sits on a
   mid-range value. Separately, 2-nonanone and LPI + 2-nonanone carry **identical** Sensory_S
   parameters (a = 7.830, b = 0.018) to three decimals, which is implausible for two independently
   fitted curves and looks like a transcription duplication. The `a` and `b` values also carry no
   uncertainty at all, so none of them can be weighted.
7. **The mass channel cannot separate nonanal from 2-nonanone.** Both are quantified at m/z 143.158
   with the same fragment at m/z 125.142. The design avoids the problem by dosing one compound per
   system, but it means **this paper can never report a mixture**, and any future use must not
   assume the two are resolvable.
8. **The stated PTR mass range is "m/z 20−25".** That cannot be right for an instrument quantifying
   at m/z 101 and 143; it is presumably a truncation of m/z 20-250 or similar. A typographic error,
   flagged so it is not propagated.
9. **The panel is ten young European women.** Non-smoking, 26 ± 2 years, saliva flow 0.145 ±
   0.1 g/min — a standard deviation of 69 % of the mean, so the panel's own salivary variability is
   large. One trained attribute ("green"), three training sessions, and the authors themselves write
   that "the variation observed in release and perception may be linked to insufficient training
   sessions" and that a carry-over effect between samples is possible. **These sensory numbers
   should never be used as an absolute perceived intensity**, only as a within-study contrast.
10. **The unflavoured lupin blank scores 51 ± 5 on the "green" scale** against 67 ± 4 for hexanal —
    the matrix supplies three quarters of the maximum flavoured percept. Figure S3 (not on disk)
    lists lupin's own descriptors as light green, grain-like, cereal, butter, fruity, barley,
    grassy, sour and lemon-like. A halo/dumping effect into the single trained attribute is likely
    and the authors say so.
11. **Water, not buffer.** The systems are MilliQ water at pH 7.0 with no buffering capacity, held
    3 h with a protein isolate at 10 g/L. **The pH is set once and never re-measured.** Every other
    binding source in the corpus buffers (Bi 2022 at 10 mM phosphate, Damodaran at 30 mM Tris,
    Leksrisompong at pH 7.0). The artificial saliva contains bicarbonate and phosphate, so the
    in-mouth pH is different again and unrecorded.
12. **The aldehyde rows carry an uncontrolled covalent share and the authors flag the chemistry
    without controlling for it.** pH 7.0 is inside the pH 6-10 window they name for Schiff-base
    formation; 3 h at 30 C is ample. There is no reducing agent, no dialysis control, no adduct
    search. **The 2-nonanone rows are the only ones in this paper free of that hazard**, which is
    why the ketone comparisons (36.84 % Imax_R suppression, 16.83 % AUC_R) are the most defensible
    numbers it contains.
13. **What this paper does NOT contain**: any binding constant in any unit; any printed binding
    percentage; any protein molar mass; any lupin site density (thiol, disulfide, amine); any
    temperature series; any pH series; any concentration series (one dose, 5 mg/L); any mixture; any
    odour threshold or OAV; any heat treatment of the protein; any measurement on pea or soy despite
    both being in the focus group; and — the authors say it themselves — any basis for
    generalisation: "Drawing conclusions about protein−aroma binding and release from exclusively
    three compounds and a simplified model system may not generalize."
14. **What to request from the authors**: (i) the nine Figure 4 binding percentages with their
    standard deviations, which is the only content in this paper that could become a registry row;
    (ii) the mucin concentration actually used, resolving the 0.01 % vs 0.216 % contradiction;
    (iii) the corrected Imax_S figure for LPI + hexanal; (iv) the intended form of eq 3;
    (v) confirmation of the Table 3 duplicate Sensory_S row for 2-nonanone; (vi) uncertainties on
    the `a` and `b` fits; (vii) a composition analysis of ProLupin 10600 beyond the 91 % / 3 %
    specification.
15. **Registry gaps.** `hexanal` and `nonanal` are keyed in `data/keys/compounds.yml`;
    **`2-nonanone` is not**, though it is already live in `parameters_matrix.py`'s
    `COMPOUND_STRUCTURE` with two shipped soy and beta-lactoglobulin constants — the product
    registry and the matrix structural registry are out of step on this compound. **Lupin protein
    isolate has no entry in `data/species/protein_matrices.yml` and no `MATRIX_LOADING` entry in
    `parameters_matrix.py`**; nor does pig gastric mucin, and mucin is not a food protein at all but
    a salivary one, so it would need a category of its own (a consumption-stage term, not a matrix
    term).
