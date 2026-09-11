# Nedvidek 1992 — EXTRACTION (alpha-dicarbonyls trapped as quinoxalines with o-phenylenediamine in situ, from glucose or xylose heated 12 h with beta-alanine or with hydrolysed wheat protein at pH 5, 6.5 and 7 — TEMPERATURE NEVER STATED; plus the identification of a new Strecker product, 5-hydroxymethyl-2-methyl-3(2H)-furanone)

### A DICARBONYL IDENTIFICATION PAPER, NOT A MELANOIDIN PAPER: it prints a fifteen-compound quinoxaline map of the sugar-fragmentation pathways, a pH 5 vs pH 7 peak-area comparison in which **1-deoxyglucosone falls 8-fold and 3-deoxyglucosone falls 39-fold while methylglyoxal is flat**, and the observation that 1,4-dideoxyosones appear only when an **alpha**-amino acid can undergo Strecker degradation — but it states no temperature for any model reaction, so nothing in it is a rate.

**Source on disk:** `data/articles/nedvidek1992.pdf` (7 pp., Z. Lebensm. Unters. Forsch. 1992,
194:222-228).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/nedvidek1992.txt`). **This is an old scanned Springer PDF and its OCR is
poor**: ligatures are mangled ("Institut fiir", "1,13 mol" appears as "1.13 mol"), Greek letters
are lost (alpha appears as "e", "~" and "a"), subscripts are corrupted, and **Table 1 and Table 2
lost structure and row labels entirely**. Pages 2, 4, 5 and 6 were therefore **rendered at
200-250 dpi (`pdftoppm`) and read from the image**; Tables 1, 2 and 3 and Scheme 2 are re-typed
below from those renders and are reliable. Where OCR remains the only source for a value, it is
marked. Schemes 1-5 are structural drawings, i.e. **figure_only** for anything numeric; Scheme 2
was rendered because the identity of compound **33** turns on it. There is no supplementary
material (1992 Springer). Repo status before this dossier: `nedvidek1992.pdf` has **no extraction
dossier** and is not cited anywhere in `src/`.

## 0. Identity

| field | value |
|---|---|
| Title | "Detection of 5-hydroxymethyl-2-methyl-3(2H)-furanone and of alpha-dicarbonyl compounds in reaction mixtures of hexoses and pentoses with different amines" |
| German title | "Nachweis von 5-Hydroxymethyl-2-methyl-3(2H)-furanon und von alpha-Dicarbonylverbindungen in Reaktionsgemischen von Hexosen und Pentosen mit verschiedenen Aminen" |
| Authors | Wolfgang Nedvidek (offprint requests), Franz Ledl — Institut für Lebensmittelchemie und Analytische Chemie; Peter Fischer — Institut für Organische Chemie; **Universität Stuttgart**, Pfaffenwaldring 55, W-7000 Stuttgart 80, Federal Republic of Germany |
| Venue | Z. Lebensm. Unters. Forsch. 1992, 194:222-228. **Received 5 September 1991.** No DOI printed (pre-DOI Springer) |
| Naming | Compounds are numbered, not named, throughout. **1-15 are quinoxalines** (the o-phenylenediamine adducts actually observed); **16-30 are the alpha-dicarbonyl compounds each quinoxaline implies**; **31-33** are synthetic intermediates and a known furanone; **34** is the new compound; **35-40** are synthetic intermediates and rejected alternatives. "Deoxyosone" means a dicarbonyl retaining the sugar's original carbon backbone; "fragmentation product" means one that does not |
| Lineage | the **Stuttgart/Munich Ledl-Severin school**. Direct predecessors: Beck, Ledl & Severin 1988 and 1989 (refs. 2, 21 — the o-phenylenediamine trapping method); Ledl & Schleicher 1990 (refs. 1, 20 — the standard review); Morita & Takagi 1986 (ref. 11 — the alkaline quinoxalines that populate column A); Severin & Seilmeier 1967 (ref. 13 — compound 33); Feather & Madson 1981 (ref. 4 — the synthetic deoxyosones); Kato 1960 (ref. 12 — the earlier bis-dinitrophenylhydrazone detection of the 3-deoxypentosone) |
| Companions on disk | `cerny2007_extraction.md`, `hofmann1998b_extraction.md` and the other Ledl-school dossiers; `martins2005_extraction.md` (whose 3-DG, 1-DG and MGO are the trunk's `TDG`, `ODG`, `MGO`); `kocadagli2016jafc_extraction.md` and `kocadagli2016foodchem_extraction.md` (the B7 furanic block's sources, which quantitate 3-DG and 3,4-DG by the same quinoxaline chemistry) |

## 1. Why it matters

This is **not** a melanoidin paper and it prints no C/N, no elemental analysis and no polymer
composition. It sits in this cluster because it maps the **dicarbonyl layer that feeds the
polymer**, and three of its compounds are trunk state variables. Concretely:

**(a) It identifies, as isolated and independently characterised quinoxalines, three of the
trunk's thirteen B1 species.** `src/kinetic_core/species.py` carries `TDG`
(3-deoxyglucosone), `ODG` (1-deoxyglucosone) and `MGO` (methylglyoxal). In this paper's numbering
those are the dicarbonyls **29** (from quinoxaline **14**), **28** (from quinoxaline **13**) and
**17** (from quinoxaline **2**). All three are found in the glucose/beta-alanine mixture, and
Table 2 gives their **relative GC peak areas at pH 5 and pH 7 side by side**. The trunk has no pH
axis on its dicarbonyl block at all; this is a directional constraint on one it might acquire.

**(b) It gives the pH direction, and it is large and not uniform.** From Table 2, going from
pH 5 to pH 7: quinoxaline 13 (**1-DG**) falls from 479 000 to 60 700 — a factor **7.9 down**;
quinoxaline 14 (**3-DG**) falls from 114 500 to 2 900 — a factor **39.5 down**; quinoxaline 2
(**methylglyoxal**) is essentially unchanged, 270 000 to 282 000 — a factor **1.04 up**; and
every small fragmentation product rises by factors of 5 to 14. The paper states the summary
itself: "**the fragmentation-product/deoxyosone ratio changes from 2:3 at pH 5 to 7:1 at pH 7**".
A model that carries 3-DG, 1-DG and MGO as a single block with one shared pH response would be
contradicted by this table.

**(c) It supplies a Strecker-degradation constraint on the amine, which is the axis this cluster
is about.** The 1,4-dideoxyosones — **24** from pentose and **30** from hexose — are found "in
much higher amounts" when an **alpha**-amino acid is present than with **beta**-alanine, and the
mechanism the paper draws (Scheme 3) is explicit: the Strecker intermediate formed from the
1-deoxyosone **reduces** it to the 1,4-dideoxyosone, releasing the Strecker aldehyde RCHO and
ammonia. **So the amine does not merely condense and leave; its decarboxylation event
re-routes the sugar's own dicarbonyl chemistry.** That is a mechanistic partner to the
decarboxylation fractions Fang 2009 and Mundt 2004 measure, and it explains why an amine's alpha
or beta constitution changes the sugar product spectrum. The trunk's `Gly` note in `species.py`
already warns that "any alpha-amino acid with the same carbon/nitrogen count substitutes without
changing the bookkeeping, **but NOT without changing the rates**" — this paper shows the effect
is not only in the rates but in **which products exist at all**.

**(d) It touches the registry twice.** Compound **33**, drawn in Scheme 2 and cited to Severin &
Seilmeier 1967, is **4-hydroxy-5-methyl-3(2H)-furanone — `norfuraneol`, which IS keyed in
`data/keys/compounds.yml`** — and this paper places it on a specific route (from the
1-deoxypentosone 22, competing with the 2,3,4-pentanetrione route). The paper's own title
compound, **34 = 5-hydroxymethyl-2-methyl-3(2H)-furanone**, is a *different* furanone and is
**not** in the registry: see Flags 8.

**(e) One thing it does NOT do, and this is the hard limit.** **No temperature is printed for any
model reaction in this paper.** The heated buffer mixtures are described only as "heated ... for
12 h"; the preparative isolation as "heated for 12 h at pH 6.5". The only temperature anywhere is
the 220 C sand-bath used to make compound 34 for detection purposes. **Nothing here can become a
rate, a barrier or a benchmark row** — see Flags 1, which is the flag that governs everything
else in this dossier.

## 2. Methods as they matter to a model

- **The trapping principle.** o-Phenylenediamine is added **to the reaction mixture itself**, not
  to a worked-up extract, and every alpha-dicarbonyl present condenses with it to a stable
  quinoxaline. The paper is candid about what this measures: the evaluation "is restricted to
  those fragmentation products which have an alpha-dicarbonyl partial structure, **or which were
  transformed into this type of compound in the course of the heating period**." So a quinoxaline
  is evidence that an alpha-dicarbonyl existed *at some point over 12 h*, not that it was present
  at any instant.
- **The four analytical mixtures.** Each was "heated in 2 ml phosphate buffer (**pH 7, 1.13 mol**
  — as printed; presumably mol/L) for 12 h", then extracted with methylene chloride and the
  residue acetylated:
  - **Mixture 1**: 60 mg (0.4 mmol) xylose, 20 mg Na2CO3·10H2O (0.07 mmol), 20 mg (0.18 mmol)
    o-phenylenediamine, **350 mg (4 mmol) beta-alanine**.
  - **Mixture 2**: the same xylose, carbonate and o-phenylenediamine, with **1 mL hydrolysed
    wheat protein** in place of the beta-alanine.
  - **Mixture 3**: 60 mg (0.33 mmol) **glucose**, same carbonate and o-phenylenediamine,
    **350 mg (4 mmol) beta-alanine**.
  - **Mixture 4**: the same glucose with **1 mL hydrolysed wheat protein**.
- **The amine.** "Hydrolyzed wheat protein (**about 4 M**) was used as a source of alpha-amino
  acids." **It is a mixture and its composition is never given** — no amino-acid profile, no
  concentration of any individual acid. Every "alpha-amino acid" result in this paper rests on
  it, except the compound-34 detection which uses pure phenylalanine. beta-alanine (4 mmol) is
  the controlled comparator because it **cannot** undergo Strecker degradation to a
  1,4-dideoxyosone in this scheme.
- **The preparative isolation** (for structure proof of quinoxalines 7, 8, 9, 10, 11): 24 g
  (148 mmol) xylose, 3 g (28 mmol) anhydrous Na2CO3, 8 g (74 mmol) o-phenylenediamine and
  **400 mL hydrolysed wheat protein**, "heated for 12 h at pH 6.5". Extracted with methylene
  chloride, then column chromatography and preparative TLC.
- **Analysis.** GC(1): Carlo-Erba 5160 Mega, FID, DB-1701 capillary 30 m x 0.32 mm x 0.25 um,
  hydrogen 40 kPa, 40 cm/s; injector and detector 270 C; 100 -> 200 C at 3 C/min, 200 -> 270 C at
  15 C/min, 30 min isothermal at 270 C. GC(2): Perkin-Elmer 8600, FID, PVMS 54 25 m x 0.32 mm x
  1 um, helium 100 kPa, 36 cm/s; 100 -> 200 C at 5 C/min, then 15 C/min to 270 C, 15 min hold.
  GC-MS: Perkin-Elmer 8420 with a Finnigan MAT Iontrap 800, **EI and CI (positive methane)**.
  NMR: Bruker 250 and 300 MHz, samples usually in CDCl3, external standard.
- **Derivatisation before GC.** Acetylation (chloroform, acetic anhydride, anhydrous sodium
  acetate, 2 h reflux) or silylation (pyridine + BSA, 1 h at room temperature). Retention times
  in the compound list are marked "(acet.)" or "(silyl)" accordingly.
- **Quantification.** "The intensities were calculated from **GC peak areas**. This
  **semi-quantitative** evaluation ..." — the authors' own word. **No internal standard is used
  for Table 2**, no response factors are applied, no replicates are mentioned, and no error is
  given. 2,3-Diphenylquinoxaline is used as an internal standard, but only in one control
  experiment (the spiking test below), not in Table 2.
- **Three control experiments on quinoxaline 11**, which matter because 11 is the paper's novel
  detection: (i) 13 mg of **32** plus 2,3-diphenylquinoxaline as internal standard were added to
  a xylose/alpha-amino acid/o-phenylenediamine mixture — "**neither a significant increase of 11,
  nor a decrease in the amount of 32 was observed**"; (ii) 11 mg of **7** was heated in 0.5 mL
  water adjusted to pH 6 with 0.1 M acetic acid for 12 h — "**compound 11 was not detectable**".
  Together these rule out 11 arising from 7 by water elimination or from 32 by oxidation during
  the run, which is what licenses the claim that the triketone **26** is genuinely formed.
- **The compound-34 detection experiment, which is a different matrix entirely.** 11 g
  phenylalanine (37 mmol) and 12 g (67 mmol) glucose were **pulverised with 13 g sea-shore sand
  and heated to 220 C for 10 min** (following Baltes & Mevissen 1988, ref. 9); volatiles were
  removed at 3 Pa and condensed in an ice trap, and the condensate examined by GC-MS against
  synthetic 34, both underivatised and silylated. **This is a dry, 220 C, 10 min, solid-diluted
  system with pure phenylalanine — it shares nothing with the buffered 12 h mixtures except the
  glucose.**
- **What is not controlled or reported.** **Temperature of every model reaction** (Flags 1); the
  composition of the hydrolysed wheat protein; how the pH 5 mixture of Table 2 was buffered (the
  methods paragraph describes only a pH 7 phosphate buffer); whether pH was re-measured; the
  atmosphere; replication; and any concentration of any product in absolute units.

## 3. Tables re-typed

### Table 1. "Quinoxalines: (A) detected in heated alkaline solutions of reducing sugars, (B) identified in glucose/beta-alanine reaction mixtures at different pH values; (C) found in xylose/alpha-amino acid reaction mixtures, (D) corresponding alpha-dicarbonyl compounds"

Read from the page-2 render. The left half of the table is drawn quinoxaline structures and the D
column is drawn dicarbonyl structures; both are written out here as formulae. **An X means the
compound was found in that system; a blank means it was not.** Column A is Morita & Takagi's
alkaline-solution result (ref. 11), i.e. **not measured in this paper**.

| # | quinoxaline (as drawn) | A | B | C | # | corresponding alpha-dicarbonyl (as drawn) | common name (mine, where standard) |
|---|---|:-:|:-:|:-:|---|---|---|
| 1 | quinoxaline (unsubstituted) | | X | X | 16 | O=CH−CH=O | **glyoxal** |
| 2 | 2-methylquinoxaline | X | X | X | 17 | O=C(CH3)−CH=O | **methylglyoxal** (pyruvaldehyde) — the trunk's `MGO` |
| 3 | 2-hydroxymethylquinoxaline | | X | X | 18 | O=C(CH2OH)−CH=O | hydroxy-2-oxopropanal (C3) |
| 4 | 2,3-dimethylquinoxaline | | X | X | 19 | O=C(CH3)−C(=O)CH3 | **2,3-butanedione** (diacetyl) — registry `2_3_butanedione` |
| 5 | 2-(2'-hydroxyethyl)quinoxaline | X | X | X | 20 | O=C(CH2−CH2OH)−CH=O | C4 fragment |
| 6 | 2-methyl-3-hydroxymethylquinoxaline | X | X | X | 21 | O=C(CH3)−C(=O)−CH2OH | 1-hydroxy-2,3-butanedione |
| 7 | 2-methyl-3-(1',2'-dihydroxyethyl)quinoxaline | X | X | X | 22 | O=C(CH3)−C(=O)−CHOH−CH2OH | **1-deoxypentosone** (1-deoxyosone, C5) |
| 8 | 2-(2',3'-dihydroxypropyl)quinoxaline | X | X | X | 23 | O=CH−C(=O)−CH2−CHOH−CH2OH | **3-deoxypentosone** (3-deoxyosone, C5) |
| 9 | 2-methyl-3-(2'-hydroxyethyl)quinoxaline | | | X | 24 | O=C(CH3)−C(=O)−CH2−CH2OH | **1,4-dideoxypentosone** |
| 10 | 2-methyl-3-vinylquinoxaline | | | X | 25 | O=C(CH3)−C(=O)−CH=CH2 | C5 enedione |
| 11 | 2-acetyl-3-methylquinoxaline | | | X | 26 | O=C(CH3)−C(=O)−C(=O)−CH3 | **2,3,4-pentanetrione** — the paper's first detection |
| 12 | 2-(1',2',3'-trihydroxypropyl)quinoxaline | | | X | 27 | O=CH−C(=O)−CHOH−CHOH−CH2OH | pentosone (the 2-oxo-aldose, C5) |
| 13 | 2-methyl-3-(1',2',3'-trihydroxypropyl)quinoxaline | | X | | 28 | O=C(CH3)−C(=O)−CHOH−CHOH−CH2OH | **1-deoxyhexosone = 1-deoxyglucosone** — the trunk's `ODG` |
| 14 | 2-(2',3',4'-trihydroxybutyl)quinoxaline | X | X | | 29 | O=CH−C(=O)−CH2−CHOH−CHOH−CH2OH | **3-deoxyhexosone = 3-deoxyglucosone** — the trunk's `TDG` |
| 15 | 2-methyl-3-(2',3'-dihydroxypropyl)quinoxaline | | X | | 30 | O=C(CH3)−C(=O)−CH2−CHOH−CH2OH | **1,4-dideoxyhexosone** — **NOT 3,4-dideoxyglucosone**, see Flags 5 |

Note the structural pattern the table encodes: **B (glucose) and C (xylose) are almost disjoint at
the top of the backbone-retaining series** — 13, 14 and 15 are hexose compounds found only in B;
9, 10, 11 and 12 are pentose compounds found only in C — while the small fragmentation products
(1-8) appear in both.

### Table 2. "Peak area of quinoxalines in glucose/beta-alanine reaction mixture of pH 5 and pH 7. Peak area obtained at"

Read from the page-4 render. **The numbers use the German thousands separator**, so "34.700"
means thirty-four thousand seven hundred. Written below with the digits regrouped, and with the
implied dicarbonyl added from Table 1 (that column is mine, not printed here).

| quinoxaline | pH 5 | pH 7 | implied alpha-dicarbonyl (mine, via Table 1) |
|---|---:|---:|---|
| 1 | 34 700 | 49 300 | glyoxal |
| 2 | 270 000 | 282 000 | **methylglyoxal (`MGO`)** |
| 3 | 3 100 | 16 200 | hydroxy-2-oxopropanal |
| 4 | 38 800 | 49 600 | 2,3-butanedione |
| 5 | 1 200 | 16 600 | C4 fragment (20) |
| 6 | 1 400 | 17 800 | 1-hydroxy-2,3-butanedione |
| **13** | **479 000** | **60 700** | **1-deoxyglucosone (`ODG`)** |
| **14** | **114 500** | **2 900** | **3-deoxyglucosone (`TDG`)** |

Quinoxalines 7, 8 and 15 are marked X in column B of Table 1 but **do not appear in Table 2**;
the paper does not say why (Flags 4). Units are bare GC-FID peak areas, uncalibrated, no internal
standard, no replicates, no error.

### Table 3. "NMR data of 5-hydroxymethyl-2-methyl-3(2H)-furanone 34, synthesized independently (see Scheme 5)"

Read from the page-6 render. Footnote as printed: "Values of delta are given relative to
tetramethylsilane; C-6 corresponds to the hydroxymethyl group, C-7 to the methyl group.
a Unresolved signal due to intermediate OH exchange rate."

| 13C | delta (ppm) | J (Hz) | value |
|---|---|---|---|
| C-3 | **205.83** | 2J(C-3, 4-H) / 2J(C-3, 2-H) / 3J(C-3, 7-H3), braced together with one arrow | **3.4** and, in parentheses, **(Σ 17.1)** |
| C-5 | **193.20** | 2J(C-5, 4-H) | 8.8 |
| | | 2J(C-5, 6-H^A,B) | 5.3 |
| | | 3J(C-5, 2-H) | 3.1 |
| C-4 | **101.80** | 1J(C-4, 4-H) | 171.1 |
| | | 3J(C-4, 2-H) | 2.2 |
| | | 3J(C-4, 6-H^A) | 2.5 |
| C-2 | **83.07** | 1J(C-2, 2-H) | 151.1 |
| | | 3J(C-2, 4-H) | 5.2 |
| | | 2J(C-2, 7-H3) | 4.6 |
| C-6 | **59.44** | 1J(C-6, 6-H^A,B) | 143.7 |
| | | 3J(C-6, 4-H) | 1.0 |
| C-7 | **16.15** | 1J(C-7, 7-H3) | 130.1 |
| | | 2J(C-7, 2-H) | 4.1 |
| **1H** | | | |
| 4-H | **5.71** | 4J(4-H, 6-H^A,B) | 1.0 |
| 2-H | **4.57** | 3J(2-H, 1-H3) | 7.2 |
| | | 5J(2-H, 6-H^A,B) | 1.0 |
| 7-H | **1.46** | 3J(7-H3, 2-H) | 7.2 |
| 6-H^A,B | — ᵃ | | |

The C-3 row is the paper's structure proof: in the fully coupled spectrum the 205.83 ppm
resonance "is split into a straightforward **sextet**", requiring coupling to **five** neighbouring
protons with virtually identical constants — 2J(4-H, C-3), 2J(2-H, C-3) and three times
3J(7-H, C-3) — which "is possible only for the proposed structure 34"; the alternative structure
**40** would show no visible coupling of its methyl protons at C-3. The Σ 17.1 is the sum over
that brace as printed.

### Every other number printed in the running text

**Reaction-mixture findings (the ones that bear on a model):**

| statement | value | where |
|---|---|---|
| **fragmentation-product / deoxyosone ratio** | **2 : 3 at pH 5; 7 : 1 at pH 7** | Results p. 226 |
| direction with pH | "distinctly higher amounts of fragmentation products are formed in the upper pH range" | Abstract; Results |
| 1,4-dideoxyosones and the amine | **24 and 30 formed in higher yields with alpha-amino acids than with beta-alanine**; "the amount of quinoxaline 9 formed **decreases significantly** when the Maillard reaction is carried out with beta-alanine instead of alpha-amino acids. The same effect is observed for the quinoxaline 15 when glucose is heated with alpha- and beta-amino acids" | Abstract; Results p. 226 — **no numbers are given for this comparison** |
| mechanism proposed | Strecker degradation of the 1-deoxyosones **22** and **28** reduces them to the 1,4-dideoxyosones **24** and **30**, releasing RCHO and NH3 | Scheme 3 |
| the small residue of 1,4-dideoxyosone without alpha-amino acids | "may be derived from the corresponding 1-deoxyosones by reaction with **reductones**" | Results p. 226 |
| 2,3,4-pentanetrione (26) | detected **for the first time**, via quinoxaline 11; proposed to form from the 1-deoxyosone **22** by cyclisation, enolisation and cleavage of the **C5-O** bond | Abstract; Results; Scheme 2 |
| the competing route from 22 | "alternative abstraction of the hydroxy group in position 2 yields the well known furanone **33**" — **read from Scheme 2 as 4-hydroxy-5-methyl-3(2H)-furanone**, i.e. `norfuraneol` (cited to Severin & Seilmeier 1967) | Results p. 226; Scheme 2 |
| the hexose analogue of the triketone | **unstable**; "spontaneous beta-diketo cleavage leads to the formation of the ester of beta-propionic acid with lactic acid" (ref. 14) | Results p. 226 |
| quinoxaline 8 | previously detected as a bis(dinitrophenyl)hydrazone (Kato 1960); **quinoxaline 7 established here for the first time** as a pentose-degradation intermediate in the presence of amines | Results p. 226 |
| **yield of compound 34** | detected in glucose/alpha-amino-acid (phenylalanine) mixtures "**in up to 0.1 % yield**" | Results p. 227 |
| pentoses in food | "In meat and meat products, pentoses are involved in the Maillard reaction" | Results p. 226 |

**Synthetic and identification data (chemistry, not kinetics):**

| compound | what is printed |
|---|---|
| **8** (2-(2',3'-dihydroxypropyl)quinoxaline) | from 100 mg (0.75 mmol) 3-deoxypentulose **23** + 100 mg (0.92 mmol) o-PD in 3 mL water, N2, room temperature, 24 h. Crystals, 50 mg, **33 %**, m.p. 96.5-99.1 C, GC(1) (acet.) t_R = 38.7 min. 13C: 154.9, 146.3, 141.4, 141.3, 130.4, 129.6, 129.3, 128.6, 71.0, 66.2, 38.2 ppm |
| **31** (2-formyl-3-methylquinoxaline) | from 430 mg (1.8 mmol) **13** + NaIO4/NaHCO3. Needles, m.p. 145.8-146.2 C, 275 mg, **74.5 %**, GC(1) t_R = 16.6 min. IR 1715, 1555, 1490, 1380, 1180, 910, 820, 770, 760 cm^-1. MS(CI) m/z 173 (M+ +1, 100 %), 144 (4) |
| **32** (2-(1'-hydroxyethyl)-3-methylquinoxaline) | from 100 mg (0.58 mmol) **31** + CH3MgBr at −78 C. Oil, 60 mg, **55 %**, GC(1) (silyl) t_R = 23.5 min. MS(CI)(silyl) m/z 262 (41), 261 (100), 246 (26), 245 (20), 216 (4) |
| **11** (2-acetyl-3-methylquinoxaline) | from 60 mg (0.32 mmol) **32** by Na2Cr2O7/H2SO4. Crystals, m.p. 87 C, 5 mg, **8.5 %**, GC(1) t_R = 19.3 min. 1H: 2.79 (s, 3H), 2.91 (s, 3H), 7.7 (m, 2H), 8.05 (m, 2H) |
| **14** (2-(2',3',4'-trihydroxybutyl)quinoxaline) | from 140 mg (0.86 mmol) **3-deoxyhexulose 29** + 115 mg (1.06 mmol) o-PD, water, N2, room temperature, 24 h. Crystals, m.p. 120 C (dec.), 80 mg, **40 %**, GC(1) (acet.) t_R = 42.8 min |
| **37** | from NaH + 5 g (31 mmol) **35** + 6.24 g (30 mmol) **36**, reflux 3 h. Oil, 1 g, **10 %** after purification, GC(2) t_R = 29.2 min |
| **38** | 1 g (3 mmol) **37**, 10 % Pd/C, H2, room temperature, pH 7, 15-60 min. Two diastereomers a and b, oil, 600 mg, **94 %**, GC(2) t_R = 19.3 and 19.9 min |
| **34** (the title compound) | 600 mg (2.8 mmol) **38** + polymer-bound pyridinium toluene-4-sulfonate, 40 C, 15-60 min. Syrupy oil, 150 mg, **42 %**. GC-MS t_R = **6.3 min**; silylated t_R = **8.23 min**. MS(EI) m/z **128 (M+, 26 %)**, 84 (31), 69 (17), 55 (100). IR 3630, 3450, 2975, 1765, 1715, 1615, 1080, 965, 800 cm^-1. UV (methanol) lambda_max **259 nm (lg eps 3.79)** |

**GC retention times and mass spectra of the quinoxalines** (all GC(1); "acet." = acetylated):
**1** 8.30 min, MS(EI) 130 (M+·, 100 %), 103 (91), 76 (69), 50 (48). **2** 11.37 min, 144 (92 %),
117 (100), 103 (8), 90 (22), 77 (42), 76 (51), 63 (11), 50 (51). **3** 26.20 min (acet.), 160
(M+−42, 100 %), 143 (6), 131 (17), 129 (20), 102 (30), 89 (3), 77 (9), 76 (22), 75 (19), 63 (7).
**4** 14.90 min, 158 (61 %), 143 (3), 117 (100), 102 (5), 90 (18), 77 (37), 76 (36), 63 (10),
50 (40). **5** 30.17 min (acet.), 217 (M+ +1, 4 %), 173 (30), 156 (100), 145 (24), 144 (20),
129 (31), 117 (7), 102 (33), 89 (9), 77 (13), 76 (33), 75 (17), 63 (9), 50 (39). **6** 29.08 min
(acet.), 216 (8 %), 174 (100), 156 (43), 143 (40), 117 (12), 116 (11), 102 (30), 89 (19), 77 (22),
76 (35), 75 (29), 63 (15), 50 (38). **7** 37.30 min (acet.), MS(CI) 289 (M+ +1, 100 %), 229 (81),
187 (31), 169 (9). **8** 38.70 min (acet.), MS(CI) 289 (100 %), 229 (4), 169 (69). **9** 32.50 min
(acet.), MS(CI) 231 (M+ +1, 67 %), 199 (11), 172 (42), 171 (100). **10** 17.16 min, MS(CI) 199
(M+ +29, 22 %), 172 (43), 171 (M+ +1, 100). **11** 19.27 min, MS(CI) 187 (M+ +1, 100 %), 158 (6).
**12** 40.90 min (acet.), MS(EI) 347 (52 %), 305 (6), 287 (100), 227 (46), 185 (12). **13**
40.30 min (acet.), MS(CI) 361 (M+ +1, 100 %), 301 (19), 241 (7). **14** 42.80 min (acet.), MS(CI)
361 (100 %), 301 (3), 241 (6), 181 (17). **15** 39.00 min (acet.), MS(CI) 303 (M+ +1, 100 %),
183 (57).

**Schemes 1-5 are structural drawings and carry no numbers; they are figure_only.**

### Arithmetic on the printed numbers (all mine)

**1. The pH 5 -> pH 7 fold changes from Table 2.** Ratio pH7/pH5, in the order the table prints
them: quinoxaline 1, **1.42x up**; 2, **1.04x up**; 3, **5.23x up**; 4, **1.28x up**; 5,
**13.8x up**; 6, **12.7x up**; **13, 7.89x DOWN**; **14, 39.5x DOWN**. The pattern is sharp: the
two backbone-retaining deoxyhexosones collapse, everything smaller rises, and methylglyoxal —
which is both a fragment and a very stable one — does not move.

**2. Does the paper's own 2:3 -> 7:1 statement reproduce from Table 2? (mine.)** Reading
"deoxyosone" as quinoxalines **13 + 14** and "fragmentation product" as **1 + 2 + 3 + 4 + 5 + 6**:
at pH 5, fragments = 349 200 and deoxyosones = 593 500, ratio **0.59 : 1**, i.e. close to the
printed **2 : 3 = 0.67**. At pH 7, fragments = 431 500 and deoxyosones = 63 600, ratio
**6.8 : 1**, against the printed **7 : 1**. **The claim reproduces from the table to within
rounding on both sides.** That is a real internal consistency check and it also tells us the
authors' fragment/deoxyosone partition is exactly the 1-6 / 13-14 split, which the paper never
states.

**3. Where the change comes from (mine).** The total peak area barely moves: 942 700 at pH 5
against 495 100 at pH 7, a factor 1.9 down. But **within** that, the deoxyosone block falls 9.3x
while the fragment block rises 1.24x. So the pH effect is **overwhelmingly a loss of the
backbone-retaining dicarbonyls**, not a gain of fragments. Any model reading this as "pH 7 makes
more fragments" would have the emphasis wrong.

**4. 3-DG versus 1-DG (mine).** At pH 5 the 1-DG quinoxaline (13, 479 000) is **4.2x** the 3-DG
quinoxaline (14, 114 500). At pH 7 it is **20.9x**. Both are uncalibrated peak areas of
*acetylated* derivatives of different molecular weight and different FID response, so **the ratio
is not a concentration ratio** (Flags 3) — but the *change* in the ratio with pH, a factor 5,
is a within-study quantity in which the response factors cancel.

**5. The compound-34 yield in context (mine).** "Up to 0.1 % yield" from 67 mmol glucose and
37 mmol phenylalanine implies at most ~37 umol of 34 if the yield is on the amino acid, or
~67 umol if on the sugar. **The paper does not say which reactant the percentage is on**, so
neither number can be quoted; it is recorded only to show the compound is a minor product.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Of everything in this paper, exactly
**three** compounds are keyed: `2_3_butanedione` (this paper's dicarbonyl **19**, via quinoxaline
**4**), `norfuraneol` (compound **33**, drawn in Scheme 2 — 4-hydroxy-5-methyl-3(2H)-furanone),
and `phenylacetaldehyde` (**not** detected here, but the Strecker aldehyde of the phenylalanine
used in the compound-34 experiment; the paper never looks for it). Glyoxal, methylglyoxal,
1-deoxyglucosone, 3-deoxyglucosone, the 1,4-dideoxyosones, 2,3,4-pentanetrione, compound **34**
and every quinoxaline are **absent** from the registry. The trunk's `TDG`, `ODG`, `MGO`, `DDG`
and `FRAG_C` are network-local names, not registry ids.

**Governing condition on every row: NO TEMPERATURE IS STATED.** Rows from Table 2 share: glucose
+ beta-alanine + o-phenylenediamine, 12 h, aqueous phosphate buffer, pH 5 or pH 7, acetylated,
GC-FID peak area, no internal standard, no replicates, semi-quantitative by the authors' own
word.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **1-deoxyglucosone quinoxaline (13), pH response** | 479 000 (pH 5) -> 60 700 (pH 7) | GC-FID peak area (arbitrary) | glucose + beta-alanine, 12 h, **T unstated** | Table 2 p. 225 | **within_study_ratio** |
| **3-deoxyglucosone quinoxaline (14), pH response** | 114 500 (pH 5) -> 2 900 (pH 7) | same | same | Table 2 | **within_study_ratio** |
| **methylglyoxal quinoxaline (2), pH response** | 270 000 (pH 5) -> 282 000 (pH 7) | same | same | Table 2 | **within_study_ratio** |
| glyoxal quinoxaline (1) | 34 700 -> 49 300 | same | same | Table 2 | within_study_ratio |
| hydroxy-2-oxopropanal quinoxaline (3) | 3 100 -> 16 200 | same | same | Table 2 | within_study_ratio |
| 2,3-butanedione quinoxaline (4) | 38 800 -> 49 600 | same | same | Table 2 | within_study_ratio |
| C4-fragment quinoxaline (5) | 1 200 -> 16 600 | same | same | Table 2 | within_study_ratio |
| 1-hydroxy-2,3-butanedione quinoxaline (6) | 1 400 -> 17 800 | same | same | Table 2 | within_study_ratio |
| **fold change, 1-DG, pH 5 -> 7** | **7.89x down** | — | same | derived from Table 2 (mine) | within_study_ratio |
| **fold change, 3-DG, pH 5 -> 7** | **39.5x down** | — | same | derived (mine) | within_study_ratio |
| fold change, MGO, pH 5 -> 7 | 1.04x up | — | same | derived (mine) | within_study_ratio |
| **fragmentation-product : deoxyosone ratio** | **2 : 3 at pH 5; 7 : 1 at pH 7** | — | same | Results p. 226 (printed); reproduces to 0.59 and 6.8 from Table 2 (mine) | **within_study_ratio** |
| 1-DG : 3-DG peak-area ratio | 4.2 (pH 5), 20.9 (pH 7) | — | same | derived (mine) | within_study_ratio (**not a concentration ratio** — Flags 3) |
| 1,4-dideoxyosone yield, alpha- vs beta-amino acid | "higher yields", "decreases significantly" — **no number** | — | glucose or xylose, hydrolysed wheat protein vs beta-alanine, 12 h, T unstated | Abstract; Results p. 226 | **level_only** (a direction with no magnitude) |
| compounds present in glucose/beta-alanine (column B) | quinoxalines 1, 2, 3, 4, 5, 6, 7, 8, 13, 14, 15 | — | as above | Table 1 | level_only (presence/absence) |
| compounds present in xylose/alpha-amino acid (column C) | quinoxalines 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 | — | xylose + hydrolysed wheat protein, 12 h, pH 6.5, T unstated | Table 1 | level_only |
| **2,3,4-pentanetrione (26)** | detected for the first time, via quinoxaline 11 | — | xylose + alpha-amino acids | Abstract; Table 1; Results | level_only (a first identification) |
| control: 11 is not an artefact of 7 | 11 mg of 7 heated 12 h at pH 6 -> "**compound 11 was not detectable**" | — | water, pH 6, 12 h, T unstated | Experimental p. 225 | level_only (a measured null) |
| control: 11 is not an artefact of 32 | 13 mg of 32 spiked in -> "neither a significant increase of 11, nor a decrease in the amount of 32" | — | full reaction mixture, internal standard 2,3-diphenylquinoxaline | Experimental p. 225 | level_only (a measured null) |
| **5-hydroxymethyl-2-methyl-3(2H)-furanone (34)** | detected, **up to 0.1 % yield** (basis unstated) | % | **glucose 67 mmol + phenylalanine 37 mmol, pulverised with 13 g sand, 220 C for 10 min**, volatiles trapped at 3 Pa | Results p. 227 | level_only (**a different matrix from every other result here**) |
| 34, mass spectrum | m/z 128 (M+, 26 %), 84 (31), 69 (17), 55 (100), EI | — | synthetic reference | Experimental p. 225 | level_only (identification) |
| 34, UV | lambda_max 259 nm, lg eps 3.79, in methanol | nm | synthetic reference | Experimental p. 225 | **measured level** (an extinction coefficient, but for a compound with no registry id and no role in any lane) |
| 34, 13C shifts | 205.83, 193.20, 101.80, 83.07, 59.44, 16.15 ppm | ppm vs TMS | CDCl3 | Table 3 | level_only (identification) |
| 34, 1H shifts | 5.71, 4.57, 1.46 ppm; 6-H unresolved | ppm vs TMS | CDCl3 | Table 3 | level_only |
| **compound 33 = norfuraneol** on the 1-deoxypentosone route | drawn as 4-hydroxy-5-methyl-3(2H)-furanone, competing with the pentanetrione route from 22 | — | pentose system | Scheme 2 (read from the render); ref. 13 | level_only (**a topology claim, not a quantity**) |
| retention times and mass spectra of quinoxalines 1-15 | see section 3 | min, m/z | GC(1), acetylated or free | Experimental p. 224-225 | level_only (method transfer) |
| synthetic yields (8, 11, 14, 31, 32, 34, 37, 38) | 33, 8.5, 40, 74.5, 55, 42, 10, 94 | % | preparative organic synthesis | Experimental | level_only (**not Maillard yields**) |
| all schemes | — | — | — | Schemes 1-5 | **figure_only** |

### How this bears on the trunk

**(a) It is not a source of constants, and the reason is Flag 1.** With no temperature, no
concentration in absolute units, no internal standard and no replicates, nothing here can be a
rate, a barrier or a benchmark row. Every usable line above is a **within-study ratio or a
presence/absence**.

**(b) Its one genuinely valuable contribution is a pH direction on three trunk species that the
trunk currently has no pH axis for.** `TDG`, `ODG` and `MGO` all live on the B1 trunk. This
paper says that between pH 5 and pH 7, in a glucose/amine pot, **3-DG collapses hardest, 1-DG
collapses substantially, and MGO does not move**. If a pH axis is ever added to the trunk, this
is the shape it has to reproduce, and the rank order (3-DG most pH-sensitive) is the testable
part.

**(c) It is a mechanistic constraint on what a Strecker event does to the sugar.** The
alpha- versus beta-amino-acid comparison shows the amine's Strecker chemistry feeds back into the
dicarbonyl pool by reducing 1-deoxyosones to 1,4-dideoxyosones. The trunk has no such edge and
this paper cannot size one (no numbers for that comparison), but it is the reason the trunk's own
`Gly` note about amine specificity is right.

**(d) It confirms an identification the B7 block depends on, and warns about a near-collision.**
The trunk's `DDG` is **3,4-dideoxyglucosone**, semi-quantitated against 3-DG in both Kocadagli
papers. Nedvidek's **30** is **1,4-dideoxyhexosone** — a different compound with a different
skeleton, a different quinoxaline (15, t_R 39.00 min, MS(CI) 303) and a different origin (Strecker
reduction of 1-DG, not dehydration of 3-DG). **They must never be matched to each other by name
similarity.** See Flags 5.

**(e) It places `norfuraneol` on a route.** Scheme 2 shows compound 33 (= norfuraneol) and the
2,3,4-pentanetrione 26 as **competing fates of the same 1-deoxypentosone intermediate 22**. The
repository keys `norfuraneol` but the trunk's B7 furanic block does not carry it; if it ever
does, this is a topology anchor — though a pentose one, and the trunk is a hexose model.

## 5. Flags

1. **NO TEMPERATURE IS PRINTED FOR ANY MODEL REACTION.** The analytical mixtures are "heated in
   2 ml phosphate buffer (pH 7, 1.13 mol) for 12 h" and the preparative run is "heated for 12 h at
   pH 6.5" — that is the whole of it. The only temperature anywhere in the paper is the **220 C**
   sand bath used to make compound 34 for detection, which is a different experiment in a
   different matrix. **This single omission is what stops every number in the paper from becoming
   a kinetic quantity.** It is the first thing to request from the authors, and until it is
   answered, Table 2's pH comparison is a shape and not a measurement at a stated condition.
2. **The pH 5 mixture is never described.** The Experimental gives one buffer, "**phosphate
   buffer (pH 7, 1.13 mol)**", for all four analytical mixtures, and the preparative run is at
   pH 6.5. Table 2 then compares **pH 5 and pH 7**. How the pH 5 mixture was buffered, whether
   the sodium carbonate was omitted or replaced, and whether pH was re-measured after 12 h are
   all unstated. The "1.13 mol" itself is printed without a volume basis.
3. **Table 2 is uncalibrated GC-FID peak areas of acetylated derivatives, and the authors call it
   semi-quantitative.** No internal standard, no response factors, no replicates, no error bars.
   Quinoxalines 13 and 14 are C6-derived and acetylated to different degrees than quinoxalines 1
   and 2; their FID responses are not equal. **Peak-area ratios ACROSS rows are not concentration
   ratios.** Ratios DOWN a row (the same compound at two pH values) are the only comparisons in
   which the response factor cancels — and those are exactly the pH fold-changes this dossier
   carries.
4. **Three compounds marked present in Table 1 column B are missing from Table 2.** Quinoxalines
   **7**, **8** and **15** carry an X in the glucose/beta-alanine column of Table 1 but have no
   row in Table 2, with no explanation. Since **15** is the 1,4-dideoxyhexosone quinoxaline —
   the compound the paper's central alpha-versus-beta claim is about — its absence from the only
   quantitative table in the paper is a real gap.
5. **1,4-dideoxyhexosone (30) is NOT 3,4-dideoxyglucosone.** The trunk's `DDG` species is
   3,4-dideoxyglucosone (3,4-DG), the rate-determining intermediate on the Kocadagli HMF limb.
   This paper's compound 30 is the **1,4**-dideoxy isomer, formed by Strecker reduction of the
   1-deoxyosone, with a methyl ketone at C1. Different skeleton, different origin, different
   quinoxaline. The names are one digit apart and the collision would be silent.
6. **The amine in every alpha-amino-acid experiment is an uncharacterised hydrolysate.**
   "Hydrolyzed wheat protein (about 4 M)" is the source, with no amino-acid profile and no
   individual concentration. So "alpha-amino acid" here means "a mixture of them at unknown
   composition", and the alpha-versus-beta comparison changes not only the amine's constitution
   but its identity, its concentration basis and its counter-ions all at once. Only the compound-34
   detection uses a single, named amino acid (phenylalanine), and that is in the 220 C sand
   experiment.
7. **The alpha-versus-beta result — the paper's headline mechanistic claim — carries no
   numbers at all.** "Higher yields", "decreases significantly", "the same effect". There is no
   table, no figure and no peak area for it anywhere. It is a direction, and it must be carried
   as `level_only`.
8. **Compound 34 is a new compound with no registry id and no established role.** The paper's own
   closing sentence is "Further experiments are necessary to obtain more information about the
   properties of 34." Its detection is in a 220 C sand-dilution experiment at "up to 0.1 % yield"
   on an unstated basis. **Do not key it or model it on this evidence.** Note that it is easily
   confused with two compounds the registry does key: `norfuraneol` (4-hydroxy-5-methyl-3(2H)-
   furanone, this paper's **33**) and `hdmf` (4-hydroxy-2,5-dimethyl-3(2H)-furanone, which this
   paper mentions only as the literature reference compound for its C,H coupling comparison,
   ref. 17). **All three are 3(2H)-furanones and all three are different.**
9. **Column A of Table 1 is not this paper's measurement.** It is Morita & Takagi's result in
   heated **alkaline** solutions (ref. 11), and the paper says outright that "because of the high
   pH value, these results are not representative of either food or biological systems". Never
   quote column A as evidence from Nedvidek.
10. **The quinoxaline trap measures a history, not a state.** The authors say the evaluation
    covers dicarbonyls "which were transformed into this type of compound **in the course of the
    heating period**". Over 12 h with o-phenylenediamine present from t = 0, a quinoxaline is a
    cumulative trap. It cannot be compared with an instantaneous concentration from a
    multiresponse fit, and Table 2 is not a set of concentrations at 12 h.
11. **OCR quality.** The text layer of this scan is bad enough that Tables 1 and 2 were
    unreadable from it and Table 3's C-3 row was corrupted. Everything in section 3 comes from
    page renders at 200-250 dpi and should be trusted over any future re-extraction of the text
    layer. Two OCR artefacts worth naming: alpha is variously rendered "e", "~", "a" and
    "R-"; and the German thousands separator in Table 2 makes "479.000" look like a decimal.
12. **What this paper does not contain**: any temperature; any rate constant; any activation
    energy; any time course (one time point, 12 h); any concentration in absolute units; any
    internal standard for the quantitative table; any replicate or error bar; any melanoidin, any
    elemental analysis and any C/N; any browning or absorbance measurement; any water-activity
    point; any single named alpha-amino acid in the buffered mixtures; any supplementary material.
13. **What to request from the authors** (Ledl died in 2006; Fischer and Nedvidek may be
    reachable, or the data may sit in Nedvidek's or Noll's Stuttgart theses — ref. 19 is
    "Noll P (1992) Thesis, University of Stuttgart (**in preparation**)"): (i) **the reaction
    temperature**, without which nothing here is usable; (ii) how the pH 5 mixture was buffered;
    (iii) peak areas for quinoxalines 7, 8 and 15, especially 15; (iv) the amino-acid composition
    of the hydrolysed wheat protein; (v) numbers behind the alpha-versus-beta comparison;
    (vi) the basis of the "0.1 % yield" for compound 34.
14. **Registry gaps against `data/keys/compounds.yml`**: only `2_3_butanedione` and `norfuraneol`
    of this paper's compounds are keyed. **Glyoxal, methylglyoxal, 1-deoxyglucosone,
    3-deoxyglucosone, the 1- and 3-deoxypentosones, the 1,4-dideoxyosones, 2,3,4-pentanetrione and
    5-hydroxymethyl-2-methyl-3(2H)-furanone are all absent.** Methylglyoxal's absence is the
    notable one: it is a B1 trunk species (`MGO`), it appears in Table 2 as the single most
    abundant quinoxaline at both pH values, and it has no registry id — the same structural gap
    the other dossiers in this cluster record, namely that the registry is a product/marker list
    and carries no Maillard intermediates.
