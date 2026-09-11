# Meynier & Mottram 1995 — EXTRACTION (glycine, lysine, cysteine or methionine at 100 mM with ribose at 70 mM in 0.2 M pyrophosphate buffer, held at pH 4.5 / 5.0 / 5.5 / 6.0 / 6.5, sealed ampoules at 140 C for 1 h; ~40 volatiles quantified against an internal standard in µg per 10 mg ribose, triplicate)

### NOT A BINDING PAPER AND NOT A THRESHOLD PAPER: this is the Reading pH-factorial that gives the corpus its sharpest constant-pH shape constraints — 2-methyl-3-furanthiol collapsing more than 152-fold and furfural more than 1 500-fold between pH 4.5 and 6.5 — and it shares only a first author with the `meynier2004` covalent-binding dossier that the matrix layer's brackets actually rest on.

**Source on disk:** `data/articles/meynier1995.pdf` (6 pp., Food Chemistry **52** (1995) 361-366).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/meynier1995.txt`). **The text layer of all four tables is unreliable** —
it interleaves the compound-name and LRI columns in Table 1 (putting furfural's 151.0 against
furylmethanol's label), drops the pH 6.5 column of several rows, and renders "tr" as `E4`, `Ot;`,
`El`, `Z3` and `;:;`. **Every one of the four tables was therefore re-read from 260 dpi page
renders of journal pages 363 and 364, and the transcriptions in section 3 are from those images,
not from the text layer.** The threshold sentence on p. 365 was likewise rendered at 300 dpi and
read from the image (Flags 6). **This paper has no figures at all** — four tables and prose. There
is no supplementary material. Repo status before this dossier: Meynier 1995 is used at second hand
by `k3_final_parameter_inventory.md` (constraint rows **B2.7**, **B2.8** and **B2.10**, and the
role table at §D.3, which assigns it **HOLD-OUT** — "directional only, and the pH axis again;
excellent as a shape test, useless as a level"), and its verdict is quoted there from
`z0_index.md`; it is named in `weerawatanakorn2015_extraction.md` as "on disk `meynier1995.pdf`,
**no dossier**". It is cited nowhere in `src/`.

## 0. Identity

| field | value |
|---|---|
| Title | "The effect of pH on the formation of volatile compounds in meat-related model systems" |
| Authors | **Anne Meynier** (corresponding) and **Donald S. Mottram** — University of Reading, Department of Food Science and Technology, Whiteknights, Reading RG6 2AP, UK. **Footnote: "Present address: INRA LEIMA, BP527, 44026 Nantes, Cedex 03, France."** |
| Venue | **Food Chemistry 52 (1995) 361-366.** Received 16 December 1993; revised version received and accepted 23 March 1994. Elsevier Science Limited, printed in Great Britain, `0308-8146/95/$9.50` |
| DOI | none printed on the scan |
| Subject | **the Maillard reaction**: how a **small** pH change (0.5 units at a time, across 4.5-6.5) changes the volatile profile of an amino-acid + ribose model system held at constant pH by a pyrophosphate buffer |
| Systems | four, run separately: **glycine**, **L-lysine**, **L-cysteine**, **L-methionine**, each with **D(-)-ribose** |
| Quantity reported | "**Approximate quantities (µg/10 mg ribose)**", by peak area against an octadecane internal standard, **with no response factors determined** — the authors call them "only approximate concentrations" |
| Structure of the paper | Tables 1-4 (one per amino acid), no figures, four discussion sections (colour and overall aroma; pyrazines and pyridines; furan derivatives; sulphur compounds) |

### How this paper relates to `meynier2004_extraction.md` (which is on disk and IS cited by the binding brackets)

| | **meynier1995** (this paper) | **meynier2004** (dossier on disk) |
|---|---|---|
| Full citation | Meynier & **Mottram**, *Food Chemistry* **52** (1995) 361-366 | Meynier, Rampon, Dalgalarrondo & **Genot**, *International Dairy Journal* **14** (2004) 681-690, doi 10.1016/j.idairyj.2004.01.003 |
| Institution | **University of Reading**, Food Science and Technology (with a Nantes "present address" footnote already attached to Meynier) | **INRA, LEIMA, Nantes** — the address the 1995 footnote points at |
| Subject | **Maillard volatile formation vs pH** | **covalent adduction of hexanal and t-2-hexenal to whey protein and sodium caseinate** |
| Chemistry | ribose + amino acid, 140 C, 1 h, sealed ampoule | aldehyde + protein, **20 C, 48 h**, aqueous, pH 6.7 |
| Instruments | SDE + GC-FID + GC-MS | front-face fluorescence, UV, PAGE, amino-acid analysis after hydrolysis |
| Output | ~40 semi-quantitative volatile levels across a 5-point pH ladder | lysine and histidine **loss**, from which `meynier2004_extraction.md` §8 derives second-order rate constants |
| Where it reaches this repository | `k3_final_parameter_inventory.md` §B.2 rows B2.7 / B2.8 / B2.10, **HOLD-OUT**; nothing in `src/` | **`src/kinetic_core/matrix_sites.py`, `BINDING_CLASSES`** — the source of the `unsaturated_aldehyde_amine` bracket (**5.3e-5 to 7.9e-5 M^-1 s^-1** at 20 C) and the upper end of the `saturated_aldehyde_amine` bracket (**≤ 2.5e-5 M^-1 s^-1**) |

**The relation is one shared first author and nothing else.** Anne Meynier did the 1995 work at
Reading with Mottram on Maillard chemistry and the 2004 work at INRA Nantes with Genot on
aldehyde-protein adducts; the 1995 paper's own footnote records the move. **A third Meynier is also
on disk** — `Meynier2002_extraction.md`, which is in fact Meynier, Garillon, Lethuaut & Genot,
*Lait* **83** (2003) 223-235, the air/skim-milk partition paper, and it is the one
`parameters_matrix.py` cites for `MATRIX_LOADING["skim_milk"]`, for five `REVERSIBLE_BINDING` rows
and for `unsat_penalty_dairy_headspace`. **So three different Meynier papers reach this repository
through three different doors, and only this one has nothing to do with the matrix layer.** The
conflation risk is real and worth stating plainly: **nothing in `meynier1995` is a binding constant,
a partition coefficient, a rate, or an odour threshold measured here.** Its numbers are Maillard
product levels.

One genuine chemical link between the 1995 and 2004 papers is worth keeping, because it runs through
the same code: `k3` row **B2.8** pairs "Zheng + Meynier" for the finding that **hydrogen-sulfide
availability rises with pH while the thiols fall with pH**. That is this paper's Table 3. The
nucleophile pools it is about — cysteine's thiol, lysine's amine — are the same pools
`src/kinetic_core/matrix_sites.py` charges as `free_thiol` and `amine` from
`data/species/protein_matrices.yml`. **The 1995 paper measures what those nucleophiles produce; the
2004 paper measures how fast they are consumed. They are two halves of the same pool's story, and
the repository currently uses only the second half.**

## 1. Why it matters

**(a) It is a constant-pH factorial, which is rare, and that is the whole point.** The paper's own
argument is that earlier pH work used unbuffered systems where "pH variation of 3 or more units may
occur during the heating", so the reported pH is not the reaction pH. Here **0.2 M pyrophosphate**
holds it: "no change of more than **0.2 pH unit** occurred after heating", verified by measuring the
pH again after the cook. Meat's own pH sits at 5.5-6.0 and moves by no more than 0.2-0.5 units in
cooking, so the design is deliberately matched to that. **This is why `k3` grades it as a shape
test**: the pH axis is real and clean even though the levels are semi-quantitative.

**(b) It is already load-bearing in the inventory's pH argument, and this dossier lets those rows be
audited.** `k3_final_parameter_inventory.md` §B.2 is the section that concludes "**NO family-level
pH term can pass**", and three of its ten constraint rows are this paper. Section 3 below re-types
all four tables from page images so those rows can be checked against the printed values — and
**one of them does not reproduce** (Flags 2).

**(c) It carries the corpus's sharpest measured pH collapse for two compounds the sulfur lane cares
about.** 2-Methyl-3-furanthiol and 2-furfurylthiol are both keyed in `data/keys/compounds.yml`
(`2_methyl_3_furanthiol`, `2_furfurylthiol`) and both appear in `matrix_oav.py`'s `_ZHOU_S2` water
threshold table (0.005 and 0.006 µg/L). This paper measures both across the same five pH values in
one pot, one panel of instruments, triplicate: **MFT 15.2 -> trace, FFT 8.5 -> 1.4** µg per 10 mg
ribose.

**(d) It offers a provenance candidate for one of the repository's uncited water thresholds.**
`_ZHOU_S2` is ingested with `provenance_flag="basis_declared_true__provenance_UNCITED"` because
Zhou 2023's SI "gives NO citation ... for any threshold in this table". This paper quotes an MFT
odour threshold from **Gasser & Grosch 1988** which, read with a unit correction that the numbers
themselves force, **brackets the repository's MFT value exactly** (Flags 6). That does not license
changing anything, but it names a paper to retrieve.

**(e) What it does NOT do.** It contains **no protein**, **no matrix**, **no binding constant**, **no
partition coefficient**, **no rate constant**, **no activation energy**, **no time course** (one
heating time only), **no absolute concentration** and **no threshold measured here**. It cannot add
a row to `MATRIX_THRESHOLDS`, cannot touch `BINDING_CLASSES`, and cannot supply a level to any fit.

## 2. Methods as they matter to a model

- **The pot, exactly.** For each system: an aqueous solution of **the amino acid at 100 mM and
  ribose at 70 mM**, made up in **0.2 M disodium dihydrogen pyrophosphate**. **Five 25 mL portions**
  into 50 mL volumetric flasks; **pH adjusted to 4.5, 5.0, 5.5, 6.0 or 6.5 with 0.2 M tetrasodium
  pyrophosphate**; volume made to **50 mL** with pyrophosphate buffer of the matching pH. So the
  amino acid and ribose are **halved to 50 mM and 35 mM** in the final mixture — the paper states
  the concentrations before the make-up and never restates them after (Flags 4).
- **Why pyrophosphate and not phosphate.** "Pyrophosphate buffer was used instead of phosphate
  buffer because of its better buffering capacity in the chosen range of pHs (4.5-6.5)." The
  authors acknowledge that phosphate **catalyses** the Maillard reaction (Potman & van Wijk 1989)
  and argue that "provided that the overall phosphate concentration remains constant, this should
  not prevent the comparison of reactions carried out at various pHs". **So every level in this
  paper is a pyrophosphate-catalysed level, and only the ratios across the pH axis are protected by
  that argument.**
- **The cook.** "The reaction mixtures were heated in **sealed glass ampoules at 140 °C under
  pressure for 1 h**" (following Whitfield 1988). **One temperature, one time. There is no kinetic
  axis in this paper at all.** Cooled to room temperature, then **the pH was measured again**.
- **The pH really did hold.** "The pyrophosphate buffer proved to be an efficient buffer in the pH
  range 4.5-6.5 since **no change of more than 0.2 pH unit occurred after heating**." This is the
  single most valuable methodological statement in the paper and it is what separates it from Shu
  1988, Zhou 2023 and every other unbuffered pH study in the corpus.
- **Volatile isolation.** **2 mL** of each reaction mixture diluted into **50 mL** of glass-distilled
  water; **Likens-Nickerson simultaneous distillation-extraction for 2 h with 20 mL diethyl ether**;
  **internal standard 100 µL of octadecane in n-hexane at 50 ng/µL (= 5 µg) added to the ether
  before extraction**; dried over anhydrous sodium sulphate; concentrated by fractional distillation
  over a Vigreux column to **~100 µL**. **SDE at atmospheric pressure for two hours is itself a
  thermal treatment**, and any artefact it generates is common to all five pH values.
- **GC.** Hewlett-Packard 5890, **on-column injector**, WCOT fused silica **30 m x 0.32 mm i.d.,
  1 µm film, DB5** (J&W); helium at **1.5 mL/min**; **40 C for 5 min, then 4 C/min to 220 C, held
  15 min**. **Linear retention indices** from a C8-C20 n-alkane solution in diethyl ether.
- **GC-MS.** Carlo Erba 4200 with split-splitless injector coupled to a **Finnigan 4000**; same
  column, direct into the ion source at **250 C**; helium at **1 mL/min**; **60 C for 5 min, then
  4 C/min to 220 C, held 15 min**; **EI at 40 eV**, continuous scan, **1 s per scan, 33-400 amu**;
  INCOS 2100 data system. Identification by spectral comparison against in-house and published
  collections (Heller & Milne 1978; ten Noever de Brauw 1980), **confirmed against authentic
  standards' LRIs "whenever possible"** — the compounds where it was not possible carry footnote a,
  "tentative identification based on comparison of mass spectra with literature spectra", and there
  are **six of them** (Flags 3).
- **Quantification, and its ceiling.** "Individual components in the volatile extracts were
  quantified by comparison of GC peak area with the area of the internal standard. **Since the
  response factors of each component were not determined, this provided only approximate
  concentrations.**" That sentence is why `z0_index.md`'s verdict on this paper reads "**Absolute
  ppb: NEVER**". A flame-ionisation detector's response per gram varies by roughly a factor of two
  to five across the classes in these tables (hydrocarbon-like thiophenes against oxygenated
  furanones), so **cross-compound comparisons within a table are worth much less than
  same-compound comparisons across the pH row**.
- **Replication and precision.** "The data were obtained from **triplicate analyses**, and
  **coefficients of variation were below 20 %**, except for a few compounds present in low
  concentrations." Every printed value is a mean of 3. **No standard deviations are printed
  anywhere.**
- **Detection limits, printed as table footnotes.** **nd = not detected, limit of detection
  c. 10 ng per 10 mg ribose**; **tr = trace, < 0.1 µg per 10 mg ribose**. So `nd` means below
  0.01 µg/10 mg and `tr` means between 0.01 and 0.1 — **the two censoring levels differ by a factor
  of ten and both are censored, which matters for every ratio that ends in `tr` or `nd`.**

## 3. Tables re-typed

**Every value below was read from a 260 dpi render of journal pages 363 and 364, not from the text
layer.** The unit throughout is **µg per 10 mg ribose**. Common footnotes on all four tables: "LRI:
linear retention index."; "nd: not detected (limit of detection c. 10 ng/10 mg ribose)."; "tr: trace
(< 0.1 µg/10 mg ribose)."; and, where marked, "a Tentative identification based on comparison of
mass spectra with literature spectra."

### Table 1. Glycine + ribose, pH 4.5 to 6.5. "Each value is the mean of 3 replicates"

| Compound | LRI | 4.5 | 5.0 | 5.5 | 6.0 | 6.5 |
|---|---:|---:|---:|---:|---:|---:|
| **Pyrazines** | | | | | | |
| Methylpyrazine | 824 | nd | nd | 0.3 | 0.5 | 0.7 |
| 2,5-dimethylpyrazine | 912 | nd | nd | 0.5 | 3.0 | 4.0 |
| **Furans** | | | | | | |
| 2-Furfural | 831 | **151.0** | 10.1 | 3.6 | 0.3 | **tr** |
| 2-Furylmethanol | 855 | 1.2 | 3.0 | 2.9 | 4.3 | 3.4 |
| 4,5-Dihydro-2-methyl-3(2H)-furanone | 805 | nd | 0.4 | 0.6 | 1.5 | 1.2 |
| 4-Hydroxy-5-methyl-3(2H)-furanone | 1042 | 12.9 | 5.2 | 5.5 | 1.2 | 0.8 |

(Table 1 carries no footnote a: every compound in it was confirmed against a standard.)

### Table 2. Lysine + ribose, pH 4.5 to 6.5

| Compound | LRI | 4.5 | 5.0 | 5.5 | 6.0 | 6.5 |
|---|---:|---:|---:|---:|---:|---:|
| **Pyrazines** | | | | | | |
| Methylpyrazine | 824 | nd | 0.4 | 1.3 | 4.1 | 4.6 |
| 2,5-Dimethylpyrazine | 912 | nd | 0.1 | 0.6 | 2.8 | 5.0 |
| **Pyridines** | | | | | | |
| 2-ethyl-4-methylpyridine (a) | 1134 | 0.9 | 1.4 | 2.5 | 3.4 | 5.4 |
| A dimethylpyridine (a) | 1262 | 0.9 | 1.2 | 1.8 | 5.6 | 6.5 |
| **Furans** | | | | | | |
| 2-Furfural | 831 | 10.7 | 4.0 | **tr** | **nd** | **nd** |
| 2-Furylmethanol | 855 | 1.4 | 1.0 | 1.4 | 1.9 | 1.9 |
| 4-Hydroxy-5-methyl-3(2H)-furanone | 1042 | 8.1 | 7.1 | 7.3 | 1.4 | 1.3 |

### Table 3. Cysteine + ribose, pH 4.5 to 6.5

| Compound | LRI | 4.5 | 5.0 | 5.5 | 6.0 | 6.5 |
|---|---:|---:|---:|---:|---:|---:|
| **Pyrazines** | | | | | | |
| Methylpyrazine | 824 | nd | nd | nd | 0.9 | 3.6 |
| Dimethypyrazine [sic] | 912 | nd | nd | nd | 0.3 | 1.5 |
| **Furans** | | | | | | |
| 2-Furfural | 831 | 3.9 | 1.5 | 1.0 | 0.6 | **tr** |
| 2-Furylmethanol | 855 | nd | 1.3 | 1.5 | 0.8 | 0.9 |
| 4,5-Dihydro-2-methyl-3(2H)-furanone | 805 | nd | 0.2 | 0.6 | 0.8 | 0.7 |
| 4-Hydroxy-5-methyl-3(2H)-furanone | 1042 | 10.6 | 6.1 | 5.2 | 2.6 | 3.0 |
| **2-Methyl-3-furanthiol** | 868 | **15.2** | 6.3 | 1.7 | 1.0 | **tr** |
| **2-Furylmethanethiol** | 909 | **8.5** | 6.0 | 3.1 | 1.5 | 1.4 |
| **Thiophenes** | | | | | | |
| 2-Formylthiophene | 1000 | 2.9 | 2.7 | 2.9 | 1.9 | 3.6 |
| 2-Formyl-5-methylthiophene | 1125 | 1.5 | 2.0 | 2.4 | 2.6 | 4.0 |
| 2-Propionylthiophene | 1188 | 3.0 | 3.7 | 3.2 | 1.5 | 1.8 |
| 4,5-Dihydro-3(2H)-thiophenone | 948 | 0.4 | 0.7 | 1.6 | 1.5 | 2.8 |
| 2-Methyl-4,5-dihydro-3(2H)-thiophenone | 987 | 5.7 | 6.7 | 6.0 | 4.8 | 5.1 |
| 2-Thiophenethiol (a) | 969 | 3.8 | 3.4 | 3.9 | 1.5 | 1.7 |
| Thieno[2,3-b]thiophene | 1205 | 1.3 | 0.8 | 0.3 | nd | nd |
| **Other sulphur compounds** | | | | | | |
| 3-Mercapto-2-pentanone | 903 | 5.4 | 5.0 | 3.3 | 0.5 | **tr** |
| Thiazole | 715 | nd | nd | 2.4 | 1.9 | 2.3 |
| 2-Acetylthiazole | 1018 | nd | 0.3 | 0.5 | 1.6 | 4.4 |
| 1,2-Dithian-4-one (a) | 1169 | 1.4 | 1.0 | 1.2 | 0.9 | 2.0 |
| 3-Methyl-1,2-dithian-4-one (a) | 1312 | tr | 0.3 | 0.7 | 0.9 | 2.0 |

(The text-layer transcription of this table was the worst of the four: it lost the pH 6.5 column of
five rows and printed "52" for 5.2, ";:;" for 2.6, "0,7" for 0.7 and "Z3" for 2.3. **The values
above are from the page image.**)

### Table 4. Methionine + ribose, pH 4.5 to 6.5

| Compound | LRI | 4.5 | 5.0 | 5.5 | 6.0 | 6.5 |
|---|---:|---:|---:|---:|---:|---:|
| **Pyrazines** | | | | | | |
| Methylpyrazine | 824 | nd | nd | nd | nd | tr |
| **Furans** | | | | | | |
| 2-Furfural | 831 | **92.4** | 10.5 | 3.2 | 2.1 | **1.9** |
| 2-Furylmethanol | 855 | 3.1 | 3.4 | 3.6 | 4.0 | 4.9 |
| 4-hydroxy-5-methyl-3(2H)-furanone | 1042 | 3.1 | 2.1 | 1.1 | 1.2 | 1.8 |
| 2-Furylmethyl methyl sulphide | 1002 | 2.5 | 3.2 | 5.5 | 7.6 | 7.5 |
| **Aliphatic sulphur compounds** | | | | | | |
| Dimethyl disulphide (a) | 727 | 16.2 | 22.2 | **24.82** | **8.6** | 41.9 |
| 3-(Methylthio)propanol [sic — the text calls it 3-(methylthio)propanal, i.e. methional] | 903 | **167.2** | 98.6 | 138.2 | 109.1 | 90.0 |
| Dimethyl trisulphide (a) | 965 | 0.5 | 0.6 | 0.6 | 0.5 | 0.5 |

**`24.82` is printed exactly so** — four significant figures where every other cell in all four
tables carries at most three (Flags 5).

### Numbers printed in the running prose

| quantity | value | where | whose |
|---|---|---|---|
| **2-methyl-3-furanthiol odour threshold** | "reported to be very low, **5-10 µg/kg**" | p. 365, verified from a 300 dpi render | **Gasser & Grosch 1988**, quoted — matrix not stated here (Flags 6) |
| pH drift after heating | **≤ 0.2 pH unit** | Results, first paragraph | this paper |
| coefficients of variation | **< 20 %**, "except for a few compounds present in low concentrations" | Results | this paper |
| cysteine / methionine pKa(NH3+) | **10.28** and **9.21** | Discussion | quoted, source not named |
| earlier component counts in these systems | "over **70** components of the glycine-ribose reaction and **120** from the cysteine-ribose" | Results | Salter 1988; Farmer 1989, quoted |
| Shu & Ho 1988 comparison | pyrazines detected at **pH 7.1** but not at **5.1 or 2.2** | Discussion | Shu & Ho 1988, quoted |
| Shu 1985 comparison | 2,5-dimethyl-4-hydroxy-3(2H)-furanone more stable at **pH 5.1** than at **2.2 or 7.1** | Discussion | Shu 1985, quoted |

**Colour and aroma are described only in words** — glycine and lysine yellow-brown at pH 4.5 and
brown at pH 6.5; cysteine and methionine yellow-orange at low pH tending red-brown at high pH;
glycine/lysine aroma caramel-like at 4.5 becoming sweet, nutty and roasted; cysteine strongly
sulphurous and unpleasant at 4.5 becoming roasted-meat-like; methionine burnt and cabbage-like at
4.5 becoming cabbage and potato-like. **No colour was measured** — there is no absorbance, no L*a*b*,
no browning index anywhere in this paper, and no sensory panel.

### Arithmetic on the printed values (all mine)

**1. The pH 4.5 -> 6.5 fold change, per compound, per system.** A `tr` cell is censored at 0.1 and an
`nd` cell at 0.01, so those ratios are **lower bounds** and are written with `>`.

| compound | glycine | lysine | cysteine | methionine |
|---|---:|---:|---:|---:|
| **2-Furfural** | **> 1 510x down** | **> 1 070x down** | **> 39x down** | **48.6x down** |
| 2-Furylmethanol | 2.8x **up** | 1.4x up | > 90x up (nd -> 0.9) | 1.6x up |
| 4-Hydroxy-5-methyl-3(2H)-furanone | 16.1x down | 6.2x down | 3.5x down | 1.7x down |
| 4,5-Dihydro-2-methyl-3(2H)-furanone | > 120x up | — | > 70x up | — |
| Methylpyrazine | > 70x up | > 460x up | > 360x up | > 1x up (nd -> tr) |
| 2,5-Dimethylpyrazine | > 400x up | > 500x up | > 150x up | — |
| **2-Methyl-3-furanthiol** | — | — | **> 152x down** | — |
| **2-Furylmethanethiol** | — | — | **6.1x down** | — |
| 3-Mercapto-2-pentanone | — | — | > 54x down | — |
| 2-Acetylthiazole | — | — | > 440x up | — |
| Thieno[2,3-b]thiophene | — | — | > 130x down | — |
| 3-(Methylthio)propanal (methional) | — | — | — | 1.86x down |
| Dimethyl disulphide | — | — | — | 2.59x **up** |
| Dimethyl trisulphide | — | — | — | **1.0x — flat** |

**2. Two of the inventory's three constraint rows reproduce exactly; the third does not.**
`k3_final_parameter_inventory.md` row **B2.10** states "MFT falls **> 152x**, FFT **6.1x**, furfural
**15-49x** over pH 4.5 -> 6.5". Against the verified tables: **MFT > 152x ✓** (15.2 -> tr, censored
at 0.1) and **FFT 6.1x ✓** (8.5/1.4 = 6.07) — both exact. **But furfural's verified range across the
four systems is 48.6x (methionine, the only uncensored one) to more than 1 510x (glycine), not
15-49x.** The upper end 49 matches methionine's 48.6; the lower end 15 matches **151.0/10.1 = 15.0,
which is glycine's pH 4.5 -> 5.0 step, not its 4.5 -> 6.5 span**. **Given that the text layer of
Table 1 mis-associates the furfural row with the furylmethanol label and drops its pH 6.5 cell —
exactly the damage this dossier had to render pages to repair — the most likely explanation is that
B2.7 and B2.10 were computed from that damaged reading.** The two rows should be re-derived from the
transcription above; the direction and the "no family-level pH term can pass" conclusion are
unaffected and if anything strengthened, since the true spread across four systems (48.6x to
> 1 510x for one compound) is **wider** than the row claims.

**3. The sharpest single statement this paper supports.** In the cysteine pot, over **one pH unit**
(4.5 to 5.5), **2-methyl-3-furanthiol falls 8.9x** (15.2 -> 1.7) while the same pot's
**2-formylthiophene is unchanged** (2.9 -> 2.9) and its **2-formyl-5-methylthiophene rises 1.6x**
(1.5 -> 2.4). **Three compounds, one pot, one pH unit, three different signs.** That is a
same-experiment sign-crossing of the kind `k3` §B.2 is built on, and it is cleaner than most of the
rows there because the pH is buffered and the instrument is one instrument.

**4. The MFT / FFT ratio inverts across the range.** At pH 4.5 the ratio is 15.2/8.5 = **1.79**; at
pH 6.0 it is 1.0/1.5 = **0.67**; at pH 6.5 it is < 0.1/1.4 = **< 0.07**. **A 25-fold swing in a
ratio of two compounds that are both cysteine-derived thiols in the same pot.** Any model that
generates both from a common H2S pool must reproduce this, and a single pH factor applied to a
"thiols" class cannot.

**5. Furfural and its alcohol run opposite ways in every system.** Glycine: furfural 151.0 -> tr
while 2-furylmethanol 1.2 -> 3.4. Methionine: 92.4 -> 1.9 while 3.1 -> 4.9. Cysteine: 3.9 -> tr
while nd -> 0.9. Lysine: 10.7 -> nd while 1.4 -> 1.9. **Four systems, four inversions, no
exceptions.** The paper offers no explanation; the obvious one — reduction of the aldehyde to the
alcohol being favoured at higher pH, or the aldehyde being consumed by amines faster than the
alcohol is — is not tested here.

**6. Amino acid identity moves furfural by 39x at fixed pH.** At pH 4.5: glycine **151.0**,
methionine **92.4**, lysine **10.7**, cysteine **3.9**. **The amino acid matters as much as the pH
does**, and the paper's own explanations are chemical and specific: lysine's extra amino group
consumes furfural into coloured products, and cysteine's hydrogen sulfide consumes it into
2-furfurylthiol. Both are sinks that a network model has to carry explicitly; **neither is a pH
effect.**

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** **Ten of this paper's compounds are
keyed**: `methylpyrazine`, `2_5_dimethylpyrazine`, `furfural` (2-furfural),
`2_methyltetrahydrofuran_3_one` (= 4,5-dihydro-2-methyl-3(2H)-furanone), `norfuraneol`
(= 4-hydroxy-5-methyl-3(2H)-furanone), `2_methyl_3_furanthiol`, `2_furfurylthiol`
(= 2-furylmethanethiol), `dimethyl_disulfide`, `dimethyl_trisulfide` and `methional`. **Absent:**
2-furylmethanol (furfuryl alcohol), 2-ethyl-4-methylpyridine, the unidentified dimethylpyridine, all
seven thiophenes, thiazole, 2-acetylthiazole, 3-mercapto-2-pentanone (note that the registry's
`mercapto_2_propanone` is a **different** compound), both dithianones, and 2-furylmethyl methyl
sulphide. **The two dithianones and the two pyridines are tentative identifications and should not
be keyed on this paper's evidence.** Every row below shares: **50 mM amino acid + 35 mM ribose after
make-up (100 mM and 70 mM before), 0.2 M pyrophosphate buffer, pH held to within 0.2 units,
140 C for 1 h in a sealed glass ampoule under pressure, SDE into diethyl ether with an octadecane
internal standard, GC-FID quantification with NO response factors determined, mean of 3 replicates,
CV < 20 %.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| every cell of Tables 1-4 (~40 compounds x 5 pH values x 4 systems) | see section 3 | **µg per 10 mg ribose** | 140 C / 1 h, pH 4.5-6.5 | Tables 1-4, pp. 363-364 (read from page images) | **level_only** — semi-quantitative by the authors' own statement; **never an absolute concentration** |
| 2-methyl-3-furanthiol, pH 4.5 -> 6.5, cysteine pot | **> 152** | x fall (15.2 -> tr, censored at 0.1) | as above | Table 3 | **within_study_ratio** — reproduces `k3` B2.10 exactly |
| 2-furfurylthiol, pH 4.5 -> 6.5, cysteine pot | **6.1** | x fall (8.5 -> 1.4) | as above | Table 3 | **within_study_ratio** — reproduces `k3` B2.10 exactly |
| 2-furfural, pH 4.5 -> 6.5 | **48.6** (methionine, uncensored); **> 39** (cysteine); **> 1 070** (lysine); **> 1 510** (glycine) | x fall | as above | Tables 1-4 | **within_study_ratio** — **does NOT reproduce `k3` B2.7/B2.10's "15-49x"** (Flags 2) |
| norfuraneol, pH 4.5 -> 6.5 | 16.1 / 6.2 / 3.5 / 1.7 | x fall, glycine / lysine / cysteine / methionine | as above | Tables 1-4 | within_study_ratio |
| methylpyrazine and 2,5-dimethylpyrazine, pH 4.5 -> 6.5 | all rise, from `nd`; **> 70x to > 500x** | x | as above | Tables 1-4 | within_study_ratio (**every ratio censored at the low end**) |
| methional, pH 4.5 -> 6.5 | **1.86** | x fall (167.2 -> 90.0), **non-monotone** | methionine pot | Table 4 | within_study_ratio |
| dimethyl disulphide, pH 4.5 -> 6.5 | **2.59** | x rise (16.2 -> 41.9), **strongly non-monotone**, with 8.6 at pH 6.0 | methionine pot | Table 4 | within_study_ratio (**the 8.6 outlier, Flags 5**) |
| dimethyl trisulphide, pH 4.5 -> 6.5 | **1.0** | x — flat at 0.5-0.6 across all five pH values | methionine pot | Table 4 | within_study_ratio — **a measured null on the pH axis** |
| MFT / FFT ratio | 1.79 (pH 4.5) -> 0.67 (pH 6.0) -> < 0.07 (pH 6.5) | — | cysteine pot | derived (mine) from Table 3 | within_study_ratio |
| furfural at fixed pH 4.5, across amino acids | 151.0 / 92.4 / 10.7 / 3.9 | µg/10 mg ribose, glycine / methionine / lysine / cysteine | 140 C / 1 h | Tables 1-4 | within_study_ratio (**38.7x by amino acid alone**, mine) |
| pH stability of the buffer | **≤ 0.2** | pH units drift over the cook | 0.2 M pyrophosphate, 140 C / 1 h | Results | level_only — **the design feature that makes the ratios usable** |
| replicate precision | **< 20 %** | CV | triplicate | Results | level_only |
| detection limit | **10** | ng per 10 mg ribose (`nd`) | GC-FID after SDE | table footnotes | level_only |
| trace ceiling | **0.1** | µg per 10 mg ribose (`tr`) | as above | table footnotes | level_only |
| 2-methyl-3-furanthiol odour threshold | **5-10** | **µg/kg as printed** — see Flags 6 | matrix not stated in this paper | p. 365, quoting **Gasser & Grosch 1988** | **threshold** — second-hand, matrix unstated, and its unit is very probably a slip (Flags 6). **Do not tabulate from this paper.** |
| absolute concentration of anything | — | — | — | — | **none exists**; the authors state response factors were not determined |
| any rate, barrier, binding constant, partition coefficient or time course | — | — | — | — | **none exists in this paper** |

### What can and cannot be done with these

**(a) They are shapes, and `k3` has already classed them correctly.** `k3_final_parameter_inventory.md`
§D.3 assigns Meynier 1995 **HOLD-OUT**, "directional only, and the pH axis again; excellent as a
shape test, useless as a level". **This dossier confirms that judgement from the primary tables and
strengthens it**: the FID response-factor gap means even the cross-compound comparisons inside one
table are soft, so what survives is the **per-compound pH profile**, which is exactly what a shape
test uses.

**(b) The one thing this dossier changes is that two inventory rows should be re-derived.** See
Flags 2. No repository value moves, because nothing from this paper is in `src/`.

**(c) It is not a matrix paper and must not be filed as one.** There is no protein, no fat, no
emulsion, no gel and no food. The "meat-related" of the title means *meat-relevant chemistry in
water*, not a meat matrix. The genuine meat-matrix paper in this batch is
`brewer1995_extraction.md`, and it is a different laboratory entirely.

**(d) It cannot be pooled with `meynier2004` or `Meynier2002` for anything.** Different chemistry,
different institution, different decade, different instruments. The only thing the three share is
Anne Meynier's name.

## 5. Flags

1. **The text layer of all four tables is untrustworthy and this dossier's transcriptions come from
   page images.** Specific failures observed: Table 1's compound labels and LRI column are
   interleaved so that **furfural's 151.0 appears against 2-furylmethanol's label**; the pH 6.5
   column is dropped from five rows of Table 3 and two of Table 1; and "tr" is variously rendered
   `E4`, `Ot;`, `El`, `Z3`, while "2.6" becomes `;:;`, "5.2" becomes `52` and "0.7" becomes `0,7`.
   **Anyone re-extracting this paper from `pdftotext` output will get several numbers wrong, and the
   wrong ones are in the furfural row that the inventory quotes.** Journal pages 363 and 364 were
   rendered at 260 dpi and every cell in section 3 was read from those images.
2. **`k3_final_parameter_inventory.md` rows B2.7 and B2.10 state furfural's pH collapse as "15-49x";
   the verified tables give 48.6x to more than 1 510x.** The MFT (> 152x) and FFT (6.1x) figures in
   the same row reproduce exactly, so the error is confined to furfural. **"15" is reproducible as
   151.0/10.1, the glycine pot's pH 4.5 -> 5.0 step**, which is consistent with the row having been
   computed from the damaged text layer described in Flags 1. **Recommend re-deriving B2.7 and B2.10
   from section 3.** The direction of the constraint and the section's conclusion do not change; the
   magnitude widens.
3. **Six compounds are tentative identifications** (footnote a): **2-ethyl-4-methylpyridine**, the
   unnamed **"A dimethylpyridine"**, **2-thiophenethiol**, **1,2-dithian-4-one**,
   **3-methyl-1,2-dithian-4-one** and **dimethyl disulphide / dimethyl trisulphide** in Table 4.
   They rest on mass-spectral comparison with literature spectra and **not** on an authentic
   standard's retention index. The dimethyl sulphides in particular are well-known compounds whose
   tentative status here is surprising and unexplained. **Do not key any of these six on this
   paper.**
4. **The stated concentrations are pre-make-up and the paper never restates them.** "an aqueous
   solution containing the amino acid (100 mM) and ribose (70 mM) was prepared in disodium
   pyrophosphate (0.2 M). Five portions of 25 ml were placed in 50 ml volumetric flasks ... The
   final volumes were then adjusted to 50 ml with pyrophosphate buffer of the appropriate pH."
   **The reacting mixture is therefore 50 mM amino acid and 35 mM ribose**, unless the 25 mL portions
   were made up from a double-strength stock, which is not stated. Note also that the yields are
   normalised **per 10 mg of ribose**, so the normalisation absorbs this ambiguity for the ratios but
   not for any attempt to reconstruct an absolute molar yield.
5. **Two anomalies in Table 4's dimethyl disulphide row.** (i) The pH 5.5 cell prints **24.82** —
   four significant figures where every other cell in all four tables has at most three; **read as
   printed and flagged, not corrected**. (ii) The row runs 16.2, 22.2, 24.82, **8.6**, 41.9: the
   pH 6.0 value is **less than half** its neighbours on both sides, in a row the abstract and the
   discussion both describe as a simple increase with pH ("as the pH increased ... an increase in the
   disulphide was observed"). **The abstract's claim is true of the endpoints and false of the
   series.** A single aberrant triplicate mean is the likely cause, but there are no standard
   deviations printed to test it.
6. **The one threshold in this paper is second-hand, its matrix is not stated, and its unit is very
   probably wrong by 1 000 — which matters because the corrected value would match a number the
   repository already carries.** Printed on p. 365 and verified from a 300 dpi render: 2-methyl-3-
   furanthiol's "odour threshold value has been reported to be very low, **5-10 µg/kg** (Gasser &
   Grosch, 1988)". **5-10 µg/kg is not "very low"** — it is a thousand times above the value
   universally quoted for this compound, and `matrix_oav.py`'s `_ZHOU_S2` carries MFT at
   **0.005 µg/L**. **Read as 5-10 ng/kg, the sentence becomes 0.005-0.010 µg/kg, whose lower bound is
   exactly the repository's uncited value.** Two readings are possible and this dossier asserts
   neither: either Meynier mis-transcribed Gasser & Grosch's ng/kg as µg/kg, or the repository's
   Zhou-derived 0.005 is 1 000x low. **Since `_ZHOU_S2` is flagged
   `basis_declared_true__provenance_UNCITED` precisely because Zhou prints no source, retrieving
   Gasser & Grosch 1988 (Z. Lebensm. Unters. Forsch. 186, 489-494, "Identification of volatile
   flavour compounds with high aroma values from cooked beef") would test both readings at once and
   might supply the missing citation.** Note the matrix caveat that would then apply: that paper is
   about **cooked beef**, so its threshold may not be an aqueous one, and a beef threshold is never
   transferred to another matrix. **Nothing should be tabulated from this paper's sentence.**
7. **The table label and the prose disagree on methional's identity.** Table 4 prints
   "**3-(Methylthio)propanol**"; the Discussion writes "**3-(methylthio)propanal (methional, the
   Strecker aldehyde from methionine)**". Methional is the **propanal**. The table's "propanol" is a
   typographical error, and the LRI (903 on DB5) is consistent with the aldehyde. **Key it as
   `methional`.** Separately, note that **3-mercapto-2-pentanone in Table 3 carries the same LRI,
   903** — the two are in different tables and different systems, so this is probably a coincidence
   rather than a mix-up, but it is worth knowing before either LRI is used to identify anything.
8. **Pyrophosphate is a Maillard catalyst and every level here is a catalysed level.** The authors
   say so, citing Potman & van Wijk 1989, and defend only the *comparison* across pH. **No level
   from this paper is comparable to a level from an unbuffered or phosphate-buffered study**, and
   the corpus contains both. This is a second, independent reason — beyond the missing response
   factors — why only ratios survive.
9. **One temperature, one time, no kinetics.** 140 C for 1 h. There is no time course, no second
   temperature, no rate and no barrier. Anything time- or temperature-resolved that appears to come
   from this paper has come from somewhere else.
10. **Every "rise with pH" ratio in this paper is censored at its low end.** Pyrazines, pyridines,
    thiazoles and the dihydrofuranone all start at `nd` or `tr` at pH 4.5, so their fold-increases
    are lower bounds set by a **detection limit that differs by 10x between the two censoring
    codes** (10 ng vs 0.1 µg per 10 mg ribose). The ">" signs in section 3 are load-bearing.
11. **No colour measurement, no sensory panel, no standard deviations.** The browning is described in
    words only ("yellow-brown", "brown", "red-brown"); the aroma likewise ("caramel-like", "sweet,
    nutty and roasted"). The CV < 20 % statement is the only dispersion information in the paper and
    it is a single global claim, not a per-value one.
12. **What this paper does not contain**: any protein, matrix, emulsion or food; any binding
    constant, partition coefficient or rate constant; any absolute concentration; any time or
    temperature axis; any threshold measured here; any water-activity variation; any figure; any
    supplementary material.
13. **What to request**: (i) the response factors, or a re-quantification against authentic
    standards, which is the only thing that would turn ~200 semi-quantitative cells into levels;
    (ii) standard deviations for the triplicates, especially for Table 4's dimethyl disulphide row;
    (iii) Gasser & Grosch 1988, to settle Flags 6.
14. **Registry gaps against `data/keys/compounds.yml`**: ten of the paper's compounds are keyed and
    listed at the head of section 4. The most consequential absences for a sulfur-lane shape test are
    the **seven thiophenes** (2-formylthiophene, 2-formyl-5-methylthiophene, 2-propionylthiophene,
    4,5-dihydro-3(2H)-thiophenone, 2-methyl-4,5-dihydro-3(2H)-thiophenone, 2-thiophenethiol,
    thieno[2,3-b]thiophene), which are the compounds whose pH profiles **cross** those of the thiols
    in the same pot and are therefore the ones that make the sign-crossing argument. **Note that
    `matrix_oav.py`'s `_ZHOU_S2` already uses the bare keys `thiophene`, `2_4_dimethylthiazole` and
    `2_thiophenecarboxaldehyde` — none of which is a registry id either**, so the thiophene family is
    a gap on both sides of the repository at once.
