# Tai & Ho 1997 (JAFC) — EXTRACTION (L-cysteine, L-cysteinesulfinic acid or L-cysteic acid + D-glucose, equimolar, water, pH 6.0 initial, sealed 120 mL bomb, 160 °C oven, 2 h; ONE end-point table of solvent-extracted volatiles in "ppm" against a tridecane internal standard)

### NOT A CYSTEINE-vs-CYSTINE PAPER. The oxidised arms are the two IRREVERSIBLE oxidation states (sulfinic and sulfonic acid), not the disulfide. Reduced cysteine + glucose gives 39 volatiles, 46 % of the total by mass sulfur-containing (thiophenes, thiazoles, cyclic polysulfides); both oxidised forms give **no sulfur volatile at all** and switch the products to pyrazines and HMF. No H2S, no free cysteine, no thiols of the MFT / FFT class, no disulfides, no time course, no replicates. Read 2026-10-31.

**Source on disk:** `data/articles/tai1997.pdf` (4 pp., born-digital ACS PDF, clean text layer;
publisher download stamp in the margin). Prose and Table 1 read from `pdftotext -layout`; **Table 1
was re-checked cell by cell against a 110 dpi render of p. 3588** (every value sits in the column the
text layer gives it). There are no figures. No Supporting Information exists for this article.

**Provenance codes:** **[M]** measured and printed · **[V]** verbatim quotation (short) ·
**[C]** cited by the authors from another paper · **[Z]** derived by this reading, never printed ·
**[Q]** qualitative statement. Nothing in this paper is digitised or fitted.

| field | value |
|---|---|
| Title | "Influence of Cysteine Oxidation on Thermal Formation of Maillard Aromas" |
| Authors | Chao-Ying Tai, **Chi-Tang Ho** (corresponding, ho@aesop.rutgers.edu) |
| Affiliation | Department of Food Science, Cook College, New Jersey Agricultural Experiment Station, Rutgers, The State University of New Jersey, New Brunswick NJ |
| Venue | J. Agric. Food Chem. **1997, 45 (9), 3586-3589**. Received 10 Apr 1997, revised and accepted 16 Jun 1997; ASAP abstract 1 Aug 1997 |
| DOI as printed | **10.1021/jf970297t** — the PDF prints no "DOI" line (1997 ACS style); it prints the manuscript code **S0021-8561(97)00297-5** and the article id **JF970297T** (p. 3589), and the publisher stamp names `jf970297t.pdf`. ✔ matches the brief |
| Companion | `tai1998_extraction.md` (same lab, same vessel / oven / 160 °C / 2 h / extraction, glutathione and glutathionesulfonic acid × pH 3 / 6 / 8) |
| Same lab, earlier | `zhang1991b_extraction.md` (glucose + Cys or GSH, 180 °C, 1 h, mg/mol) |

## 1. What this paper answers for the model

**The brief asked for** measured yields or time courses of thiols (MFT, FFT, 2-mercapto-3-pentanone),
thiophenes, thiazoles, H2S and disulfides (bis(2-methyl-3-furyl) disulfide etc.), their dependence on
cysteine vs **cystine** and on pH, with conditions and quantitation. What the paper actually holds:

| asked for | in this paper? | what there is |
|---|---|---|
| cysteine vs **cystine** | **NO** | the oxidised arms are **cysteinesulfinic acid** (RSO2H) and **cysteic acid** (RSO3H). Cystine is named only in the introduction [Q]. The brief's description of the paper is wrong on this point |
| MFT, FFT, 2-mercapto-3-pentanone | **NO** | none reported in any arm. The only thiol-bearing small molecules are **2-mercapto-3-butanone 1.54 ppm** and **1,1-ethanedithiol 0.66 ppm** (reduced Cys arm) [M] |
| thiophenes, thiazoles | **yes, one end point, one arm** | 12 thiophenes (sum 20.41 ppm), 3 thiazoles (sum 3.00 ppm) in Cys + glucose [M]/[Z] |
| H2S | **NO** | invoked as the explanation ("the lack of hydrogen sulfide"), never measured [Q] |
| disulfides | **NO** open-chain disulfide; cyclic polysulfides only | 7 cyclic polysulfides (trithiolane, trithiane, two tetrathianes, pentathiepane, dithiolane, ethanedithiol), sum 31.70 ppm [M]/[Z] |
| pH dependence | **NO** | one initial pH (6.0), never re-measured |
| time course | **NO** | one end point (2 h at oven setpoint 160 °C) |
| free cysteine / residual precursor | **NO** | — |
| quantitation | **semi-quantitative** | FID area against **one** internal standard (tridecane), response factors not stated (implicitly 1); unit "ppm" with **no stated basis**; no replicates, no SD |

**For the cysteine-excess question (free Cys 2-40× too high in the engine's fed-TTCA pots;
`kinetic_core_dpo_partition_prereg.md` secs. 2, 7):** the paper measures neither cysteine nor any
oxidation product formed in situ, so it **cannot discriminate a missing sink from oxidation**. It
supports two narrower statements:

1. **Volatiles are not the missing sink.** Even reading "ppm" as µg per g of the 50 g pot (the most
   generous plausible basis), all 23 sulfur volatiles together carry **≈ 44 µmol S = 0.18 % of the
   24 mmol cysteine sulfur** (arithmetic 3) [Z]. With an unknown response factor this could be off by
   a factor of a few, not by the ~two orders of magnitude needed to matter.
2. **If a pot's cysteine were oxidised past the disulfide (sulfinic / sulfonic), its sulfur
   chemistry would stop:** zero sulfur volatiles from either oxidised form, the authors' reading
   being that these forms release no H2S [M]/[Q]. That is a **property of a fed oxidised precursor**;
   whether sulfinic / sulfonic acid forms in a sealed Maillard pot is not tested here. The
   oxidation channel relevant to the engine (thiol → disulfide, `OX` / `ch_cys_ox`, B11) is the one
   this paper does **not** study.

**Fit rows: none. Hold-outs: none worth two-reader keying** (§4). The one transferable item is a
directional claim (oxidised sulfur donor → no sulfur volatiles, more pyrazines / HMF), which the
engine cannot be scored on because it carries no sulfinic or sulfonic species.

## 2. Methods as they matter to a model [M]

- **Charge.** "L-Cysteine or L-cysteic acid (**0.024 mol**) ... with an equal amount of α-D-glucose,
  was dissolved in **50 mL** of distilled water" → **0.48 mol/L each** [Z]. For cysteinesulfinic acid +
  glucose, "**0.002 mol** of each reactant was used" → **0.04 mol/L each** [Z], **12× lower** than the
  other two arms; the volume is not restated (assumed 50 mL).
- **pH.** "adjusted to **pH 6.0**" — reagent not named, **unbuffered**, final pH never measured. Cysteic
  acid (a strong sulfonic acid) needs about one extra equivalent of base to reach pH 6, so the three
  arms differ in ionic strength and counter-ion load [Z].
- **Vessel and heating.** "sealed in a **120 mL bomb reactor**", "heated at **160 °C in an oven for 2 h**",
  then ice-water quench. Headspace ≈ 70 mL of air [Z]; vessel material and come-up time not stated
  (an oven-heated steel bomb lags its setpoint by tens of minutes; the 2 h is not 2 h at 160 °C).
- **Volatile isolation.** Whole reaction mixture + **1 mL internal standard (tridecane, 100 ppm)**;
  liquid-liquid extraction with CH2Cl2 **5 × 50 mL**; Na2SO4; Kuderna-Danish to 5 mL; N2 to 1 mL.
- **GC.** Varian 3400, FID, **DB-1 60 m × 0.25 mm × 1 µm**, split 25 : 1, 40 → 280 °C at 2 °C/min,
  He 1 mL/min. **GC-MS**: Varian 3400 + Finnigan MAT 8230 (EI 70 eV). **Identification by mass
  spectra against the NBS library and earlier papers only — no retention indices** (contrast Tai
  1998, which adds RI).
- **Sensory.** Odour descriptions by two flavorists [Q].
- **Replication.** None stated; no error of any kind printed.

## 3. Table 1 re-typed (verbatim). "Volatile Compounds Generated from Thermal Reaction of Glucose with Cysteine, Cysteinesulfinic Acid, or Cysteic Acid" — amount, ppm [M]

Columns: **CG** cysteine / glucose; **CSFAG** cysteinesulfinic acid / glucose; **CAG** cysteic acid /
glucose. A blank cell is printed blank (not detected / not reported).

| class | compound | CG | CSFAG | CAG |
|---|---|---|---|---|
| thiophenes | thiophene | 1.50 | | |
| | 2-methylthiophene | 1.91 | | |
| | 2-ethylthiophene | 0.61 | | |
| | 2-acetylthiophene | 7.78 | | |
| | 2,3-dimethylthiophene | 0.54 | | |
| | 2,5-dimethylthiophene | 0.74 | | |
| | thieno(3,2-b)thiophene | 2.56 | | |
| | tetrahydrothiophen-3-one | 1.06 | | |
| | 2-methyltetrahydrothiophen-3-one | 1.52 | | |
| | 5-methyltetrahydrothiophen-3-one | 0.33 | | |
| | 5-methyl-2-thiophenecarboxylic acid | 1.45 | | |
| | 2-butanoylthiophene | 0.41 | | |
| thiazoles | 2-methylthiazole | 0.97 | | |
| | 2,5-dimethylthiazole | 0.11 | | |
| | 2,4,5-trimethylthiazole | 1.92 | | |
| polysulfides | 1,1-ethanedithiol | 0.66 | | |
| | 2-methyl-1,3-dithiolane | 0.85 | | |
| | 3,5-dimethyl-1,2,4-trithiolane | 11.42 | | |
| | 3-methyl-1,2,4-trithiane | 9.94 | | |
| | 3,6-dimethyl-1,2,4,5-tetrathiane | 2.60 | | |
| | 4,6-dimethyl-1,2,3,5-tetrathiane | 2.61 | | |
| | 4,7-dimethyl-1,2,3,5,6-pentathiepane | 3.62 | | |
| pyrazines | pyrazine | | 5.90 | 4.49 |
| | methylpyrazine | | 9.84 | 15.69 |
| | 2,3-dimethylpyrazine | | 1.57 | 2.12 |
| | 2,5-dimethylpyrazine | | 22.90 | 30.59 |
| | ethylpyrazine | | 4.12 | 1.40 |
| | vinylpyrazine | | 0.36 | |
| | trimethylpyrazine | | 11.01 | |
| | 2-ethyl-5-methylpyrazine | | 6.95 | 8.88 |
| | 2-ethyl-6-methylpyrazine | | 4.80 | 2.00 |
| | 2-vinyl-5-methylpyrazine | | | 0.67 |
| | 2-vinyl-6-methylpyrazine | | | 0.22 |
| | 2-ethyl-3,5-dimethylpyrazine | | 3.04 | 1.43 |
| | 3-ethyl-2,5-dimethylpyrazine | | 27.38 | 2.04 |
| | 2-methyl-5-isopropylpyrazine | | | 0.62 |
| | 2,3-diethyl-5-methylpyrazine | | 2.38 | |
| | 3,5-diethyl-2-methylpyrazine | | 1.48 | |
| | 2,5-dimethyl-3-propylpyrazine | | | 1.61 |
| | quinoxaline | | | 0.49 |
| pyrroles | 2-acetylpyrrole | | 0.77 | 5.04 |
| | 5-methyl-2-formylpyrrole | | | 0.61 |
| furans | 2-ethylfuran | 0.66 | | |
| | 2,5-dimethylfuran | | 0.17 | |
| | 2,5-dimethyl-3(2H)-furanone | | 2.59 | |
| | 2-acetylfuran | 7.18 | 1.36 | 3.26 |
| | 2-furfural | | 0.82 | 0.98 |
| | 5-methyl-2-furfural | | | 8.63 |
| | 5-(hydroxymethyl)-2-furfural | | 18.66 | 36.26 |
| | 2-furanmethanol | 12.44 | 1.34 | 4.41 |
| | 1-(2-furyl)-1,2-propanedione | 0.61 | | |
| | 2,3-dihydro-5-methylfuran | 26.86 | | |
| carbonyls | 2-butanone | 6.28 | | 0.85 |
| | 1-hydroxy-2-propanone | 1.88 | 28.02 | 5.08 |
| | acetoin | 0.77 | 8.43 | 26.32 |
| | 1-methoxy-2-propanone | | | 0.12 |
| | 3-methyl-2-butanone | 0.94 | | |
| | 3-hexanone | 1.00 | | |
| | 2-mercapto-3-butanone | 1.54 | | |
| | 2-methylcyclopentanone | 0.31 | | |
| | 3-methyl-2-cyclopenten-1-one | 1.27 | | 0.22 |
| | 2-methyl-3-heptanone | 1.84 | | |
| | 3-methyl-2,4-pentanedione | 1.03 | | |
| | cyclotene | 2.41 | | 0.05 |
| | 2,3-dihydro-3,5-dihydroxy-6-methyl-4H-pyran-4-one | | 5.49 | 1.71 |
| miscellaneous | acetic acid | 0.67 | 4.11 | 0.74 |

Footnotes as printed: a cysteine/glucose; b cysteinesulfinic acid/glucose; c cysteic acid/glucose.

### Arithmetic [Z] (`scratchpad tai_scripts/tai_arith.py`, run in the maillard env; not shipped)

1. **Compound counts.** CG **39** rows, CSFAG **24**, CAG **29**. The text prints 24 and 29 (match) but
   "a total of **38** compounds" for CG (Flag 2).
2. **Column totals.** CG 122.80, CSFAG 173.49, CAG 166.53 ppm. CG by class: polysulfides 31.70,
   thiophenes 20.41, thiazoles 3.00, 2-mercapto-3-butanone 1.54 → **sulfur volatiles 56.65 ppm = 46 %
   of CG**; furans 47.75 (of which 2,3-dihydro-5-methylfuran 26.86 and 2-furanmethanol 12.44);
   other carbonyls 17.73.
3. **Sulfur recovery, assuming ppm = µg per g of the ≈ 50 g pot.** Σ (ppm × 50 g / MW × n_S) over the
   23 sulfur compounds = **43.9 µmol S** against 24 000 µmol cysteine = **0.18 %**. Largest single
   items: 3,5-dimethyl-1,2,4-trithiolane 3.75 µmol (0.016 mol % of Cys), 2-acetylthiophene 3.08 µmol
   (0.013 %), 2,4,5-trimethylthiazole 0.75 µmol, 2-mercapto-3-butanone 0.74 µmol (0.003 % each). The
   basis is an assumption (Flag 1); the conclusion "volatiles are < 1 % of the sulfur" survives any
   plausible basis and any FID response factor within ×5.
4. **The oxidised arms, per mole of precursor.** Pyrazines sum 101.73 (CSFAG) vs 72.25 ppm (CAG); if
   ppm is a mass fraction of the pot, CSFAG's 12× lower charge makes its pyrazine yield **≈ 17×** CAG's
   per mole, and HMF ≈ 6× [Z]. The authors compare the columns **as printed** and never normalise
   (Flag 3).
5. **1-hydroxy-2-propanone / acetoin**: 3.32 (CSFAG), 0.193 (CAG) — the paper's "approximately 3:1"
   and "1:5" [M] reproduce.
6. **Headspace oxygen budget.** 70 mL air × 0.2095 / 24.47 L/mol = **0.60 mmol O2** in the bomb;
   with 4 thiol per O2 (thiol → disulfide) that can oxidise at most **2.4 mmol = 10 %** of the 24 mmol
   cysteine. In this pot oxidation by the headspace cannot be a major cysteine sink whatever its rate.

## 4. What the repository can use

| item | value | conditions | source | evidence class | use |
|---|---|---|---|---|---|
| sulfur volatiles from Cys + Glc, by class | polysulfides 31.70, thiophenes 20.41, thiazoles 3.00, 2-mercapto-3-butanone 1.54 ppm | 0.48 M each, pH 6.0 initial, 160 °C oven, 2 h, sealed | Table 1 [M] | semi_quantitative_single_IS, end point, basis unstated | **level-only context**; not a fit row |
| sulfur recovered in volatiles | **≈ 0.18 %** of cysteine S | as above; ppm read as µg/g | arithmetic 3 [Z] | derived_assumption | bound: volatile products are not a cysteine sink |
| oxidised Cys (RSO2H, RSO3H) + Glc | **zero** sulfur volatiles; pyrazines 101.7 / 72.3 ppm; HMF 18.7 / 36.3 ppm | 0.04 / 0.48 M, same | Table 1 [M] | directional (presence / absence) | **directional claim** only; the engine has no sulfinic / sulfonic species, so nothing scores it today |
| reduced Cys + Glc makes **no pyrazines and no HMF** | absent | same | Table 1 [M] | directional | qualitative cross-check of the engine's H2S-vs-NH3 competition story; not scorable without a pyrazine lane |
| odour notes | CG "roasted, meaty, soy sauce-like"; CSFAG "strong coffee"; CAG "smoky, nutty, yeasty" | — | Results [Q] | qualitative | none |

**Two-reader keying: not recommended.** One end point per precursor, an undefined concentration basis,
precursor charges that differ 12×, an incomplete sulfur list (Flag 4), no replicate or error, and no
species the engine carries except as a lumped class. If a sulfur-class signature is ever wanted from
this lab, `tai1998_extraction.md` Table 1 is the better source (stated unit, three pH arms, RI-confirmed
identities).

## 5. Flags

1. **"ppm" has no stated basis.** Not per g of pot, per g of precursor, or relative to the internal
   standard — nothing says. 1 mL of "100 ppm" tridecane is ~100 µg only if the solvent's density is ~1
   (solvent not named). Every absolute reading of Table 1 rests on an assumption; ratios **within a
   column** are the only safe use.
2. **Compound count.** CG has 39 rows; the text says 38.
3. **Unequal charges compared as equal.** CSFAG was run at 0.002 mol (0.04 M), CG and CAG at 0.024 mol
   (0.48 M); the discussion ("pyrazines are produced dramatically in tremendous amount" for CSFAG) reads
   the ppm columns as if the charges matched.
4. **The sulfur list is incomplete by the authors' own later statement.** Tai & Ho 1998 (p. 2263):
   2,5-dimethyl-4-hydroxy-3(2H)-thiophenone (thiofuraneol) was formed in this cysteine / glucose system
   but "we did not identify this compound until after the study was published", so it is absent from
   Table 1. In the 1998 glutathione / glucose pots it is the **largest single volatile** (240-259 mg/mol
   at pH 6-8). The CG sulfur sums of §3 are therefore lower bounds.
5. **"3,5-Dimethyl-1,2,4-trithiolanes, the most predominant compounds in this study"** [V] — not so in
   CG by the table's own numbers: 2,3-dihydro-5-methylfuran (26.86) and 2-furanmethanol (12.44) exceed
   the trithiolane (11.42). Most predominant *sulfur* compound is what the table supports.
6. **"No or only trace amounts" (abstract) vs none (Table 1, text):** the table lists no sulfur
   compound in either oxidised arm; the abstract hedges. Cysteic acid *alone* (no glucose) gave a
   "slightly sulfury smell" [Q], unquantified.
7. **No final pH, no vessel material, no come-up time, no residual precursor, no replicates.**
8. **Identification by MS library match only** (no RI); several assignments (2,3-dihydro-5-methylfuran
   at 26.86 ppm as the largest CG peak; 2-butanoylthiophene) are unusual for a cysteine / glucose pot
   and unconfirmed.
9. **160 °C** is above the 90-150 °C range of the sulfur-lane pots the engine is fitted on; an
   oven-heated bomb's real thermal history is unknown.
10. **Registry gaps:** cysteic acid, cysteinesulfinic acid, 2-mercapto-3-butanone, the cyclic
    polysulfides and thiofuraneol are absent from `data/keys/compounds.yml`. Paper not in
    `papers.yml` before this reading.
11. **Literature lead for the cystine question** (the one this paper does not answer): Mulders 1973,
    Z. Lebensm.-Unters. Forsch. 152:193-201, "Volatile compounds from the non-enzymatic browning
    reaction of the cysteine/cystine-ribose system" [C] — cited here; DOI not printed and not verified
    in this reading.
