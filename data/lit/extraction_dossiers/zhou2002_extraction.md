# Zhou 2002 — EXTRACTION (2-pentylpyridine against commercial soy protein isolate, purified 7S beta-conglycinin and purified 11S glycinin at 0.5 % w/w = 5 g/L; equilibrium dialysis across a 3 500 MWCO membrane with headspace-SPME/GC-MS isotope-dilution readout; bound amount at pH 4.5/7/9, at 4/25/74 C, at three NaCl levels and under UV; Klotz constants at 25 C only)

### THE ONLY AQUEOUS PROTEIN-FLAVOUR BINDING MEASUREMENT ABOVE 60 C IN THIS BATCH — AND IT GOES THE WRONG WAY. Wave B26's record says "37 C is an in-mouth temperature, not a process one: nothing here licenses a pea binding constant at 90 or 140 C." This paper dialyses soy protein against 2-pentylpyridine **in water at 74 C** and reports the bound amount beside a 25 C control on the same protein: **SPI 0.396 against 0.159, 7S 0.469 against 0.235, 11S 0.569 against 0.375 moles per mole of protein (Table 1, p. 143)** — binding **RISES 1.5-2.5x on heating from 25 to 74 C**, against the dry-phase soy result (Aspelund: adsorption *falls* 1.2-1.8x per 10 C) and against the physical-adsorption expectation. The paper's explanation is thermal denaturation exposing sites, the same mechanism Guo 2019 and Crowther 1980 invoke. **But no binding CONSTANT is reported at 74 C** — the Klotz constants (SPI **107 ± 10**, 7S **131 ± 16**, 11S **228 ± 30 M^-1**) are all at 25 C, they are fitted on a 0.3-3.0 mM branch that is 10x above the concentration Table 1 was measured at, and **they do not reproduce Table 1's own 25 C values, falling short by 6-17x (mine, section 3)**. Separately: **2-pentylpyridine received a compound-registry id in this repository today, and this paper is the first source in the corpus that would put anything in its `seen_in` list** — it is also the paper that ties `2_pentylpyridine` to `e_e_2_4_decadienal`, the other registry id in its formation pathway.

**Source on disk:** `data/articles/zhou2002.pdf` (4 pp., J. Food Sci. 67 (1) 2002, 142-145).
Read from the `pdftotext -layout` text layer, with **Table 1 (p. 143) verified cell-by-cell against the rendered
page image** — necessary because the text layer rendered every degree sign as "8C" and every ± as "6". The
re-typed table below is from the page image. **Figure 1 (p. 144) is the only figure**: binding curves of bound
2-pentylpyridine (mole/mole protein, y-axis 0-30) against 2-pentylpyridine concentration (mM, x-axis 0-3) for
the 7S fraction, the 11S fraction and whole SPI, with error bars. **Its points are figure-only and are not
typed here**; note that the saturation levels it shows are far above the n values printed in the text, which is
one of two internal inconsistencies recorded in section 3. There is **no supplementary material**. Repo status
before this dossier: Zhou 2002 has **no extraction dossier** and is **not cited** in
`src/kinetic_core/parameters_matrix.py`, `src/kinetic_core/matrix_sites.py` or
`data/species/protein_matrices.yml`. **Note the name collisions on disk**: `zhou2000_extraction.md`,
`zhou2023_extraction.md`, `zhou2024_extraction.md`, `zhou2025_extraction.md` and `zhou2025b_extraction.md`
already exist and are different papers; this dossier is for `zhou2002.pdf`, J. Food Sci. 67(1):142-145.

## 0. Identity

| field | value |
|---|---|
| Title | "Binding Properties of 2-Pentyl Pyridine to Soy Protein as Measured by Solid Phase Microextraction" |
| Authors | **A. Zhou, W. L. Boatright, L. A. Johnson, M. Reuber**. Per the credit box on p. 145: *"Author Zhou is affiliated with Ross Products Division, Abbott Laboratories, Columbus, Ohio. Authors Boatright and Reuber are affiliated with the Animal Science Dept., University of Kentucky, Lexington, KY 40546-0215. Author Johnson is with the Center for Crops Utilization Research, Iowa State University, Ames. Direct correspondence to author Boatright (E-mail: wlboat1@pop.uky.edu)."* |
| Venue | Journal of Food Science, **Vol. 67, Nr. 1, 2002, pp. 142-145**. Section: "JFS: Food Chemistry and Toxicology". © 2002 Institute of Food Technologists |
| DOI | **NO DOI IS PRINTED IN THE PDF.** The manuscript-tracking line on p. 145 reads "MS 20001537 Submitted 10/26/00, Revised 5/22/01, Accepted 5/30/01, Received 5/30/01" and nothing DOI-shaped appears on any of the four pages. |
| Funding | USDA NRICGP Grant # KY9701687. Published with the approval of the Director of the Kentucky Agricultural Experiment Station as Journal Article Number 00-07-178 |
| Proteins (three) | **(1)** Commercial soy protein isolate **Supro 500E**, Protein Technologies International, St. Louis — **protein content not stated**. **(2)** **beta-conglycinin (7S)** and **(3) glycinin (11S)** fractions produced at the Center for Crops Utilization Research, Iowa State, by scaling up Nagano's method (Nagano 1992; Wu 1999). **7S fraction: 97.7 % protein, of which 70.4 % is beta-conglycinin. 11S fraction: 92.1 % protein, of which 90.6 % is glycinin.** All on a dry basis |
| The ligand | **2-Pentylpyridine (2-pp)**, Lancaster (Windham, NH). Deuterium-labelled 2-pentylpyridine (d5/d6) made after Boatright 1998 and used as an internal standard |
| Molar-mass conventions used | **SPI 100 000** (O'Neill & Kinsella 1987), **beta-conglycinin 160 000** (Thanh & Shibasaki 1978), **glycinin 320 000** (Iyengar & Ravestein 1981) — all three **cited, not measured**, and the SPI figure is explicitly called "hypothetical" |
| Naming | "binding properties" in Table 1 = **the moles of 2-pentylpyridine bound by one mole of protein**, i.e. the Scatchard v; "equilibrium binding constant" = the Klotz K in M^-1; "primary binding sites" = n |
| Companions on disk | `damodaran1981_extraction.md` (**the dialysis method this paper modifies, and the source of three shipped soy rows on the same 100 000 g/mol basis**), `guo2019_extraction.md` (soy preheat, aqueous, 37 C — the same two-class picture), `crowther1980_extraction.md` and `aspelund1983_extraction.md` (dry soy; **Crowther 1980 is cited by name in this paper's introduction**), `bi2022_extraction.md` (pea, Klotz/headspace) |

## 1. Why it matters

**1. It is the temperature point the registry does not have.** `REVERSIBLE_BINDING` carries 21 rows and their
measurement temperatures are 25, 30, 37 and 40 C. Wave B26's own note draws the line explicitly: 37 C is an
in-mouth temperature and "nothing here licenses a pea binding constant at 90 or 140 C". Of the five papers in
this batch, **Zhou 2002 is the only one that measures binding IN WATER at a temperature above 60 C**. Aspelund
1983 reaches 80-100 C but on dry protein by gas-solid chromatography; Crowther 1980 measures at 60-80 C but
also dry; Guo 2019 and Xu 2022 preheat at 80-100 C but measure at 37 C or below. **Zhou dialyses a soy protein
slurry against 2-pentylpyridine in a 74 C water bath and reports what is bound.** That is the point of contact
with the process temperature the model runs in.

**2. And its sign contradicts the dry-phase evidence.** From 25 C to 74 C the bound amount **rises**:
**SPI ×2.49, 7S ×2.00, 11S ×1.52 (mine, Table 1)**. Aspelund's dry soy loses 1.2-1.8x of its adsorption per
10 C rise over 90-100 C and calls that "characteristic of a physical adsorption process". **The two soy results
run in opposite directions over overlapping temperature ranges**, and the reason is not a contradiction but a
mechanism difference the corpus now has evidence for on both sides: in water, heating a globular protein
*denatures* it and creates surface, and the created surface outweighs the thermodynamic penalty; on a dry
surface there is nothing to denature and only the penalty acts. **A model that extrapolated a binding constant
above 60 C using a van 't Hoff slope from any dry-phase source would get the SIGN wrong on an aqueous soy
matrix.** That is the single most important thing this paper says to the registry, and it is a warning rather
than a parameter.

**3. But it does not close the gap, and the reason must be stated precisely.** Table 1's temperature block is a
**bound amount at a single ligand concentration (0.03 mM)**, not a binding constant. The paper's only binding
constants — 107, 131 and 228 M^-1 — are at **25 C**, and the abstract's own careful phrasing keeps the two
apart: *"More 2-pp was also bound at high temperature (74 °C) than at 25 °C, but greater binding affinity of
2-pp was observed at 4 °C than at 25 °C."* **There is no printed constant at 4 C or at 74 C anywhere in the
paper**, so the "affinity" half of that sentence has no number behind it (Flags 3). What can be carried is a
**within-study ratio of bound amounts**, which — in the linear (low-occupancy) regime the 0.03 mM dosing sits
in — is proportional to nK and therefore is a defensible relative temperature coefficient. **What cannot be
carried is an absolute aqueous soy binding constant at 74 C. It does not exist in this paper.**

**4. The compound is one this repository just registered, and this paper is its only binding source.**
`data/keys/compounds.yml` line 430 carries `id: 2_pentylpyridine` with SMILES `CCCCCc1ccccn1`, InChIKey
`HSDXVAOHEOSTFZ-UHFFFAOYSA-N`, `aliases: []`, **`seen_in: []`** and `identity_source: structure entered by
hand`. It has a chemical identity and **no observation attached to it anywhere in the corpus**. This paper
supplies: an aqueous soy binding constant (107 ± 10 M^-1 at 25 C), a per-gram form on the registry's own
100 000 g/mol soy basis (**3.21e-3 L/g, mine**), a fraction-resolved breakdown (7S and 11S), a pH dependence,
a salt dependence, a temperature dependence and a UV dependence. **It also ties `2_pentylpyridine` to
`e_e_2_4_decadienal`** — the other registry id in the pathway, which is also carried with `seen_in: []`. The
mechanism is printed on p. 142: **ammonia condenses with 2,4-decadienal to a Schiff base, ring closure gives a
dihydropyridine, oxidation gives 2-pentylpyridine** (Buttery 1977), **confirmed by the same authors with
13C-2,4-decadienal and 15N-ammonia** (Zhou & Boatright 2000). That is a lipid-oxidation-to-nitrogen-heterocycle
channel, and it is exactly the kind of cross-lane link the compound registry exists to make findable.

**5. It is a fraction-resolved measurement, which the registry has none of for soy.** The three shipped
Damodaran soy rows and the Arai-via-Damodaran row are all on whole isolate. This paper measures **the same
ligand against whole SPI, against purified 7S and against purified 11S under identical conditions**:
K = 107 / 131 / 228 M^-1 and n = 3 / 7.5 / 9. **Glycinin binds 2-pp with 2.1x the affinity of whole isolate and
1.7x that of beta-conglycinin (mine).** It also flags that this ordering is ligand-dependent — O'Neill &
Kinsella 1987 found the *opposite* order for 2-nonanone, and the paper says so (p. 144). **So a whole-isolate
constant is a composition-weighted average whose weighting differs by compound**, which is a real caution for
any attempt to transfer a soy constant to a differently-fractionated preparation.

**6. It gives the registry three environmental sensitivities it currently carries for almost nothing.**
`MatrixParameter` has `ph_of_measurement` as a first-class field, and only Leksrisompong's caseinate rows carry
a measured pH dependence (diacetyl binds at pH 7, not at pH 5.5). Zhou gives a **three-point pH series on three
proteins**: bound amount rises **1.8-2.8x from pH 4.5 to pH 9 (mine)**, monotonically, on all three, with the
chemically obvious explanation that 2-pp is a base and is protonated and released at low pH. It also gives a
**three-point NaCl series**: raising NaCl from 0.171 to 1.711 M **cuts binding to 0.22-0.54 of its value
(mine)**. **There is no ionic-strength term anywhere in `parameters_matrix.py`**, and this is a measurement of
how large one would be — a factor of 2-4.5 over a salt range a real formulation spans.

**What this paper does NOT give the repository**: any binding constant at any temperature other than 25 C; any
thermodynamic quantity (no ΔH, ΔS or ΔG anywhere — so the activation-energy category error cannot arise here);
any rate constant; any covalent-adduct evidence; any aldehyde, ketone or ester; any 2-alkenal (so **nothing
touches `ALPHA_BETA_UNSATURATION`**); any homologous series (2-pp is the only ligand, so **nothing checks
`CHAIN_LENGTH_SLOPE_PER_CH2`**); any protein content for the commercial isolate; any thiol, disulfide or amine
assay; any preheat-then-measure design (the 74 C arm heats *during* the measurement, which is a different and
in some ways better thing — and a worse one, Flags 1).

## 2. Methods as they matter to a model

- **The three proteins.** Commercial **Supro 500E** SPI (Protein Technologies International) — **no protein
  content is given for it anywhere**. **7S beta-conglycinin: 97.7 % protein, 70.4 % of the fraction is
  beta-conglycinin. 11S glycinin: 92.1 % protein, 90.6 % of the fraction is glycinin**, both dry basis, both
  made at Iowa State by a scaled-up Nagano method. **Note the 7S fraction is only 70 % pure in its target
  protein** — nearly a third of it is something else (Flags 7).
- **The pot.** **Equilibrium dialysis.** Spectra/Por tubing, **id 45 mm, molecular-weight cut-off 3 500**.
  **Equal volumes of SPI aqueous slurry (0.5 %, w/w) OUTSIDE the tubing and 2-pp methanol solution (0.03 mM)
  INSIDE.** So the protein loading in its own compartment is **5 g/L**, and the ligand starts entirely on the
  other side of the membrane and diffuses across. Shaken **14 to 20 h** to reach equilibrium. **Control: an
  identical dialysis in pure water at the same temperature.**
- **Loading and dosing, and the volume ambiguity that follows.** 5 g/L in the protein compartment; **0.03 mM
  2-pp initially in the ligand compartment**. Because the two compartments are of equal volume and the ligand
  equilibrates across, **the 2-pp concentration averaged over the whole system after equilibration is 0.015 mM**,
  and the protein averaged over the whole system is 2.5 g/L. The paper's Table 1 footnote says only *"The
  initial concentration of 2-pp in protein slurry was 0.03 mM"*, and **it never states which volume basis v is
  computed on** (Flags 4). For the binding *curves* the ligand ranged **0.3 to 3.0 mM**.
- **Measurement family: EQUILIBRIUM DIALYSIS with HEADSPACE-SPME QUANTITATION — a hybrid, and the paper
  explains exactly why.** It is a modification of Damodaran & Kinsella 1981a's dialysis and O'Keefe 1991's
  headspace method. Damodaran's route extracts the ligand with an organic solvent after dialysis; that failed
  here — *"Chloroform and methylene chloride failed to effectively extract 2-pp from protein slurries because
  proteins precipitated. Other solvents, such as ethanol, ethyl ether, and isooctane, could not efficiently
  separate the aqueous phase from the organic phase"*, and ultracentrifuging at **70 409 × g for 3 h** worked
  for SPI but not for the 7S and 11S fractions (p. 143). O'Keefe's static headspace lacked detection limit and
  purge-and-trap foamed the protein. **So the separation is dialysis and the readout is headspace SPME.**
  **For the registry's `method` field this is `equilibrium_dialysis`** — the free ligand is physically separated
  by a membrane, which is what the k2 sec. B.3 method boundary is really about — **but with two provenance
  qualifications that matter**: the quantitation is a **headspace ratio**, and **there is no thiol-blocking
  agent**. Damodaran's shipped `kg_nonanal_soy` row carries the note "DIALYSIS + 2-mercaptoethanol: this
  constant EXCLUDES the cysteine-aldehyde chemistry a headspace determination would count." **Zhou's dialysis
  has no 2-mercaptoethanol, no Tris, no azide — just water.** For a pyridine this matters less than for an
  aldehyde (2-pp has no carbonyl and is not a Michael acceptor), but it should be recorded.
- **The bound-concentration equation.** **B = T × (I − O) / I**, where B is the bound 2-pp concentration, T the
  total 2-pp concentration, **I the 2-pp headspace concentration inside the membrane tubing** and **O the
  headspace concentration outside** (p. 143). This is a **ratio of two headspace measurements**, so an absolute
  SPME calibration offset cancels — the same reason Amendment 4 lets Meynier's and Leksrisompong's suspect
  absolute static-headspace scales cancel in a within-run ratio. **But it assumes the headspace above the two
  compartments responds identically**, and the protein compartment contains 5 g/L of a surface-active protein
  while the other does not (Flags 5).
- **The SPME readout.** **50 mL aliquots** from each side into a **150 mL glass bottle**, septum-sealed; high
  concentrations diluted with nanopure water; **deuterated 2-pp (d5/d6) added as an internal standard**, its
  concentration kept **within 5x of the unlabelled 2-pp**. **Polydimethylsiloxane SPME fibre exposed to the
  headspace for 1 h**, desorbed **4.5 min at 210 C, 1 min splitless**, on a Hewlett-Packard G1800A GCD with EI
  detection (240 C). Column **DB-225, 30 m × 0.25 mm, 0.25 µm**, held 40 C for 5 min, ramped 3 C/min to 155 C,
  then 20 C/min to 210 C. Scan **m/z 10-200**, 5 min solvent delay, helium 1.0 mL/min. **Quantitation by the
  ratio of the m/z 93 ion to the m/z 99 and 98 ions of the internal standard** — an **isotope dilution assay**
  after Schieberle & Grosch 1987, chosen "to minimize any changes of 2-pp during the entire sampling process".
  **Separate standard curves for the inside and the outside of the membrane.** Triplicate per treatment.
- **The four treatment blocks, and note they are NOT under one common condition.**
  - **pH**: SPI solutions made in **1 M phosphate buffers at pH 4.5, 7.0 and 9.0**, at **25 C**. One molar
    phosphate is an extremely high ionic strength — and the paper's own salt block shows salt suppresses
    binding (Flags 2).
  - **Temperature**: the semipermeable-membrane dialysis was **held in a water bath at 4, 25 or 74 C**, with a
    pure-water control at the same temperature. **No buffer is mentioned for this block.**
  - **NaCl**: **0.171, 0.855 and 1.711 M** in the SPI slurries, at 25 C.
  - **UV**: the protein slurry exposed **during the entire dialysis** to direct radiation from a **dual
    wavelength UV lamp (365 and 254 nm)**.
  **So the pH block is in 1 M phosphate and the other three blocks are not**, and the four blocks' own 25 C
  reference values differ — SPI reads 0.141 (pH 7 in buffer), 0.159 (25 C), 0.217 (0.171 M NaCl) and 0.154
  (UV control). **Four nominally comparable controls spanning 0.141 to 0.217, a 1.54x spread (mine)**
  (Flags 2).
- **Protein quantitation.** **Bradford dye-binding**, calibrated against **bovine serum albumin**, "assuming
  that all protein constituents have equal binding affinity for Coomassie brilliant blue" — the paper states
  the assumption. **No protein was found inside the membrane tube after dialysis** (a real and useful control:
  the membrane held).
- **Statistics.** ANOVA in SAS; LSD at P < 0.05; means compared by **Tukey-Kramer HSD**. Table 1 footnote:
  means ± standard deviation, **n = 3**, letters within a row series.
- **The binding curves and the Klotz fit.** Ligand **0.3 to 3.0 mM**. *"Saturation was reached at 1.8 mM of
  2-pp for the glycinin fraction, and 2.4 mM for the beta-conglycinin fraction. No obvious saturation was
  observed for whole SPI even at 3 mM of ligand"* (p. 144). **According to the Klotz equation (1949)**, K and
  the primary binding sites were calculated: **K = 107 ± 10, 131 ± 16, 228 ± 30 M^-1 and n = 3, 7.5, 9** for
  SPI, 7S and 11S. **The n values are printed without uncertainties.** All at 25 C.

## 3. Tables re-typed

Table 1 was verified against the rendered page image. Evidence marks: `[M]` measured in this study, `[C]` cited
from elsewhere, `[F]` fitted.

### Table 1 (p. 143). "Binding properties^d of 2-pentyl pyridine^e to commercial soy protein isolates, 7S soy protein fraction and 11S soy protein fraction"

Values are **moles of 2-pentylpyridine bound per mole of protein**, ± standard deviation, n = 3.

| | SPI | 7S | 11S |
|---|---|---|---|
| **pH effects** | | | |
| pH 4.5 | 0.119 ± 0.002 a `[M]` | 0.208 ± 0.002 a `[M]` | 0.297 ± 0.003 a `[M]` |
| pH 7 | 0.141 ± 0.002 b `[M]` | 0.284 ± 0.003 b `[M]` | 0.366 ± 0.003 b `[M]` |
| pH 9 | 0.328 ± 0.003 c `[M]` | 0.401 ± 0.005 c `[M]` | 0.522 ± 0.006 c `[M]` |
| **Temperature** | | | |
| 4 °C | 0.241 ± 0.005 a `[M]` | 0.310 ± 0.003 a `[M]` | 0.417 ± 0.005 a `[M]` |
| 25 °C | 0.159 ± 0.003 b `[M]` | 0.235 ± 0.003 b `[M]` | 0.375 ± 0.003 b `[M]` |
| **74 °C** | **0.396 ± 0.005 c** `[M]` | **0.469 ± 0.006 c** `[M]` | **0.569 ± 0.005 c** `[M]` |
| **NaCl** | | | |
| 0.171 M | 0.217 ± 0.002 a `[M]` | 0.310 ± 0.002 a `[M]` | 0.385 ± 0.004 a `[M]` |
| 0.855 M | 0.105 ± 0.002 b `[M]` | 0.203 ± 0.002 b `[M]` | 0.299 ± 0.004 b `[M]` |
| 1.711 M | 0.048 ± 0.001 c `[M]` | 0.126 ± 0.002 c `[M]` | 0.206 ± 0.002 c `[M]` |
| **UV** | | | |
| Control | 0.154 ± 0.003 a `[M]` | 0.226 ± 0.003 a `[M]` | 0.353 ± 0.002 a `[M]` |
| Exposure | 0.353 ± 0.005 b `[M]` | 0.429 ± 0.005 b `[M]` | 0.501 ± 0.005 b `[M]` |

Footnotes exactly as printed: *a-c: Means (± standard deviation) within a row series with no common superscripts
differ (P < 0.05; n = 3).* *d: Binding properties were expressed as the moles of 2-pentyl pyridine that were
bound by one mole of SPI or soy protein fraction.* *e: The initial concentraion [sic] of 2-pp in protein slurry
was 0.03 mM.*

**Read "within a row series" as "within a BLOCK": the letters run down each of the four blocks within a column,
not across the three proteins.** (Every value in a block carries a distinct letter in all four blocks and all
three columns, so every difference within every block is significant at P < 0.05 — an unusually clean result
given standard deviations of 1-2 % of the value.)

### Constants printed in the running text, p. 144

| quantity | SPI | 7S (beta-conglycinin) | 11S (glycinin) | class |
|---|---|---|---|---|
| equilibrium binding constant K, **at 25 °C** | **107 ± 10 M^-1** | **131 ± 16 M^-1** | **228 ± 30 M^-1** | `[F]` (Klotz fit on the 0.3-3.0 mM curve) |
| primary binding sites n | **3** | **7.5** | **9** | `[F]` (no uncertainty printed) |
| molecular weight used | **100 000** ("hypothetical") | **160 000** | **320 000** | **`[C]`** — O'Neill & Kinsella 1987, Thanh & Shibasaki 1978, Iyengar & Ravestein 1981 |
| saturation concentration in the binding curve | none observed to 3 mM | **2.4 mM** | **1.8 mM** | `[M]` |
| protein content of the preparation | **not stated** | 97.7 % (70.4 % beta-conglycinin) | 92.1 % (90.6 % glycinin) | `[M]` |

### Numbers printed in the running text

| quantity | value | where | class |
|---|---|---|---|
| protein loading | **0.5 % (w/w)** slurry = **5 g/L** | p. 143 | `[M]` |
| initial ligand concentration in Table 1 | **0.03 mM** 2-pp, in methanol solution | p. 143 and Table 1 footnote e | `[M]` |
| ligand range for the binding curves | **0.3 to 3.0 mM** | p. 143 | `[M]` |
| dialysis membrane | Spectra/Por, **id 45 mm, MWCO 3 500** | p. 143 | `[M]` |
| equilibration time | **14 to 20 h** | p. 143 | `[M]` |
| internal standard | deuterated 2-pp (d5/d6), kept **within 5x** of the unlabelled concentration | p. 143 | `[M]` |
| SPME | PDMS fibre, headspace **1 h**, desorb **4.5 min at 210 °C** | p. 143 | `[M]` |
| quantitation ions | **m/z 93** (2-pp) against **m/z 99 and 98** (internal standard) | p. 143 | `[M]` |
| aliquot / vial | **50 mL into a 150 mL bottle** | p. 143 | `[M]` |
| ultracentrifugation that was tried and abandoned | **70 409 × g for 3 h** — precipitated SPI but not the 7S or 11S fractions | p. 143 | `[M]` (a measured null) |
| membrane integrity control | **"No protein was found inside the membrane tube after dialysis"** | p. 143 | `[M]` |
| **the odour of 2-pentylpyridine** | **"penetrating-grassy aroma when detected by GC-Olfactometry and throat-catching taste in water when evaluated by sensory panelists"** | p. 142 | `[C]` (Boatright & Crum 1997) |
| **the formation mechanism** | **ammonia condenses with 2,4-decadienal → Schiff base → ring closure → dihydropyridine → oxidation → 2-pentylpyridine** | p. 142 | `[C]` (Buttery 1977) |
| **the isotopic confirmation** | both precursors confirmed with **13C-2,4-decadienal and 15N-ammonia** in SPI processing | p. 142 | `[C]` (Zhou & Boatright 2000) |
| **amino acids that INCREASE 2-pp** | **arginine, lysine, asparagine, glutamine** | p. 142 | `[C]` (Zhou & Boatright 2000) |
| **amino acids with NO effect on 2-pp** | **aspartic acid, glutamic acid, glycine, histidine** | p. 142 | `[C]` (Zhou & Boatright 2000) |
| what raises 2-pp levels | pro-oxidants **FeCl3, CuCl2**, or **UV exposure** of the protein slurry | p. 142 | `[C]` (Boatright 1998) |
| the pH explanation | 2-pp is **basic**, so "it would be more easily released under acidic conditions than under alkaline conditions" | p. 144 | interpretation |
| the temperature explanation | more bound at 74 C "due to **thermal denaturation that enhanced soy protein ability to bind 2-pp**" | p. 144 | interpretation |
| the 4 C explanation | Damodaran & Kinsella found **increased hydrophobicity at 5 C compared to 20 C**; lower temperature "allowed protein unfolding, either increasing the strength or number of binding sites" | p. 144 | `[C]` + interpretation |
| the salt explanation | salt "could increase the vapor pressure" and "destabilizes electrostatic interactions, thus decreasing flavor binding" | p. 144 | interpretation |
| **the UV control that WAS run** | *"the overall 2-pp contents in the headspaces for samples exposed to UV light in this equilibrium study were similar to the control, it is unlikely that the increased binding observed in this study was the result of 2-pp synthesis"* | p. 144 | **`[M]` — a mass-balance check, and note it was run ONLY for the UV block (Flags 1)** |
| the UV explanation | chemical oxidation of the protein with "changes in the surface hydrophobicity of the proteins that are associated with protein oxidation" | p. 144 | interpretation (Boatright & Hettiarachchy 1995) |
| the contrary literature | **O'Neill & Kinsella 1987 found beta-conglycinin has a GREATER affinity for 2-nonanone than glycinin** — the reverse of the 2-pp order found here | p. 144 | `[C]` |
| the consistent literature | O'Keefe 1991 found the glycinin > beta-conglycinin order for "some carbonyl flavor compounds" | p. 144 | `[C]` |

**Figure-only quantities.** Every point, error bar and saturation level in Figure 1 (p. 144). The y-axis runs
0-30 mole/mole protein and the three curves reach values far above the printed n of 3, 7.5 and 9 (section 3,
internal check 2). Per house rule no figure value is typed here.

### Arithmetic on the printed constants (all mine)

**1. The temperature coefficient, as a ratio (mine).** Bound amount at one temperature over another, from
Table 1's temperature block:

| protein | 4 C / 25 C | **74 C / 25 C** | 74 C / 4 C |
|---|---:|---:|---:|
| SPI | 1.52x | **2.49x** | 1.64x |
| 7S beta-conglycinin | 1.32x | **2.00x** | 1.51x |
| 11S glycinin | 1.11x | **1.52x** | 1.37x |

**Binding is NON-MONOTONE in temperature on all three proteins: it falls from 4 to 25 C and then rises steeply
to 74 C**, with 74 C the highest of the three in every case. **The 25 C point is a minimum**, which is the one
temperature at which every other number in this paper was measured.

**Why the ratio is meaningful even though the absolute scale is not.** The Scatchard isotherm is
v = n[L]/(1 + K[L]). At [L] = 0.03 mM = 3e-5 M and K ≤ 228 M^-1, **K[L] ≤ 6.8e-3, so v ≈ nK[L] to better than
1 %** — the linear, low-occupancy regime. **In that regime v is directly proportional to nK, so a ratio of two
v values at the same [L] is a ratio of two nK values.** Provided the free ligand concentration is the same in
both arms (it is, to the extent that the same total dose was used and the bound fraction is small), **the
temperature ratios above are ratios of binding constants.** That is the strongest statement this paper supports
about temperature, and it is a relative one.

**2. The per-gram constants at 25 C, on the registry's own construction (mine).** K_g = nK / MW, exactly as
`kg_2_heptanone_soy` = n·K/100 000 was built from Damodaran:

| protein | n | K, M^-1 | nK, M^-1 | MW used | **K_g = nK/MW, L/g (mine)** |
|---|---:|---:|---:|---:|---:|
| whole SPI | 3 | 107 ± 10 | 321 | 100 000 | **3.21e-3** |
| 7S beta-conglycinin | 7.5 | 131 ± 16 | 982.5 | 160 000 | **6.14e-3** |
| 11S glycinin | 9 | 228 ± 30 | 2 052 | 320 000 | **6.41e-3** |

**The SPI molar mass, 100 000, is the same convention Damodaran used and the same one three shipped soy rows
rest on.** So `kg_2_pentylpyridine_soy` = 3.21e-3 L/g would sit on **exactly** the footing of
`kg_2_heptanone_soy` (4.40e-3 L/g) — same protein species, same molar-mass convention, same method family
(dialysis), same aqueous phase, 25 C in both. **That is the cleanest cross-paper basis match in this batch**,
and it puts 2-pentylpyridine just below 2-heptanone in per-gram soy affinity. Note also that **the 7S and 11S
per-gram values agree to 4 % with each other despite a 1.7x difference in K** — the molar-mass difference
(160 000 vs 320 000) almost exactly cancels the affinity difference, which is what one expects if the two
fractions have similar site density per gram and differ mainly in oligomer size.

**3. The internal check that FAILS: the printed K does not reproduce Table 1's own 25 C values (mine).** In the
linear regime, v = nK[L]. At [L] = 0.03 mM:

| protein | nK, M^-1 | **v predicted at 0.03 mM (mine)** | **Table 1, 25 C** | **discrepancy** |
|---|---:|---:|---:|---:|
| SPI | 321 | 0.0096 | 0.159 | **16.5x** |
| 7S | 982.5 | 0.0295 | 0.235 | **8.0x** |
| 11S | 2 052 | 0.0616 | 0.375 | **6.1x** |

**Table 1's bound amounts are 6-17x larger than the paper's own binding constants predict.** And the
discrepancy is a **lower bound**, because the free ligand concentration is at most the total: with SPI at
5 g/L and MW 100 000 the protein is 5e-5 M, so v = 0.159 corresponds to 7.95e-6 M of bound 2-pp out of a
3e-5 M total (mine) — about a quarter of it — which would push the free concentration down and the implied
constant up further. **The most likely reading, and it matches Guo 2019 exactly, is that there is a
high-affinity class at low concentration that the Klotz fit on the 0.3-3.0 mM branch does not see.** Guo found
primary constants 10-70x the secondary ones on soy for esters, from the same kind of split. **The implied
low-concentration nK from Table 1 is 5 300 (SPI), 7 833 (7S) and 12 500 M^-1 (11S), i.e. K_g of 5.30e-2,
4.90e-2 and 3.91e-2 L/g (mine, all upper-bounded on total ligand).** **These two sets of numbers are 6-17x
apart and the paper offers no reconciliation. Do not ship either without a decision on which branch is meant**
(Flags 6).

**4. The second internal check that fails: n against Figure 1's saturation (mine, qualitative).** The printed
primary site counts are 3, 7.5 and 9, while Figure 1's y-axis runs to 30 mole/mole protein and all three curves
rise well above 9 over the 0.3-3.0 mM range. **A site count of 3 cannot be the saturating capacity of a curve
that reaches into the teens.** The consistent reading is again a two-class picture in which only the primary
class is reported — the same structure Guo 2019 makes explicit. The paper does not say this. (No figure value
is quoted here; the statement rests only on the printed axis range against the printed n.)

**5. If one nevertheless carries the 25 C constant to 74 C using this paper's own ratio (mine, and heavily
caveated).** K_g(25 C) × [v(74 C)/v(25 C)]:

| protein | K_g at 25 C | **implied K_g at 74 C (mine)** |
|---|---:|---:|
| whole SPI | 3.21e-3 L/g | **8.00e-3 L/g** |
| 7S | 6.14e-3 L/g | **1.23e-2 L/g** |
| 11S | 6.41e-3 L/g | **9.73e-3 L/g** |

**This is a two-step chain and both steps are compromised**: the 25 C constant is the one that fails check 3,
and the 74 C ratio is confounded by 14-20 h of thermal exposure with an unchecked synthesis pathway (Flags 1).
**These three numbers are recorded so that a later reader does not recompute them and mistake them for the
paper's; they are not shippable.** What is shippable in spirit is the *direction and rough size*: **an aqueous
soy matrix binds roughly 1.5-2.5x MORE at 74 C than at 25 C.**

**6. pH, salt and UV as ratios (mine).**

| effect | SPI | 7S | 11S |
|---|---:|---:|---:|
| pH 9 / pH 7 | **2.33x** | 1.41x | 1.43x |
| pH 7 / pH 4.5 | 1.19x | 1.37x | 1.23x |
| **pH 9 / pH 4.5** | **2.76x** | **1.93x** | **1.76x** |
| 0.855 M / 0.171 M NaCl | **0.48x** | 0.66x | 0.78x |
| **1.711 M / 0.171 M NaCl** | **0.22x** | **0.41x** | **0.54x** |
| **UV exposure / control** | **2.29x** | **1.90x** | **1.42x** |

**Three observations.** (a) **Whole SPI is the most sensitive protein to every one of the four factors** —
largest pH swing, largest salt suppression, largest UV effect, largest temperature effect — and 11S glycinin
the least, in all four. That is a consistent pattern across sixteen independent comparisons and it says the
commercial isolate's binding is dominated by a labile, surface-dependent component that the purified fractions
have less of. (b) **The salt effect is large enough to matter to a formulation**: 1.711 M NaCl (10 % w/v)
removes 78 % of whole-SPI binding. `parameters_matrix.py` has no ionic-strength term at all. (c) **The pH
effect on SPI, 2.76x over 4.5 units, is comparable in size to the temperature effect over 70 C** — so a
`ph_of_measurement` field is doing real work for this compound, which is exactly why the registry made it
first-class.

**7. The fraction ordering (mine).** On K at 25 C: **11S / 7S = 1.74x**, **11S / SPI = 2.13x**, **7S / SPI =
1.22x**. The paper says "the glycinin fraction showed a 2-fold greater affinity for 2-pp than the
beta-conglycinin fraction"; the printed constants give **1.74x**, not 2.0x. On the Table 1 bound amounts at
25 C the ordering is the same but the ratios differ: 11S/7S = 1.60x, 11S/SPI = 2.36x (mine). **And on the
per-gram basis the ordering collapses**: 3.21e-3, 6.14e-3, 6.41e-3 L/g — 7S and 11S are indistinguishable and
whole SPI is half of both (section 3 item 2). **Which fraction "binds best" depends entirely on whether one
asks per mole or per gram, and only the per-gram question is the one the model asks.**

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** **`2_pentylpyridine` is keyed** (id `2_pentylpyridine`,
line 430; SMILES `CCCCCc1ccccn1`, InChIKey `HSDXVAOHEOSTFZ-UHFFFAOYSA-N`, `identity_source: structure entered
by hand 2026-09-01`) and its **`seen_in` list is EMPTY** — this paper is the first corpus source that would
populate it. **`e_e_2_4_decadienal` is also keyed** (line 712) and also carries `seen_in: []`; this paper is
the link between the two. **Neither is in `COMPOUND_STRUCTURE` in `parameters_matrix.py`**, which carries no
nitrogen heterocycle of any kind — adding 2-pentylpyridine would mean adding a class (`alkylpyridine`), not
just a key, exactly as Wave B26 found for the first alcohol. **2-Pentylpyridine is not on
`data/species/off_flavour_targets.yml`**, so no scored panel consumes it today and no prediction changes on
its account.

Every row below shares: **soy protein at 0.5 % w/w = 5 g/L in its compartment, dialysed across a 3 500 MWCO
membrane against 0.03 mM 2-pentylpyridine in methanol for 14-20 h, headspace SPME on a PDMS fibre with
deuterated-2-pp isotope dilution, GC/MS, triplicate.** Unless stated otherwise the medium is **water**, not a
buffer, and the temperature is **25 C**.

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **Klotz binding constant K, 2-pentylpyridine × whole soy protein isolate** | **107 ± 10** | M^-1 (per mole of protein, MW 100 000 stated) | **25 C**, water, 5 g/L, fitted on 0.3-3.0 mM | p. 144 | **`binding_constant`** |
| Klotz binding constant K, 2-pentylpyridine × 7S beta-conglycinin | **131 ± 16** | M^-1 (MW 160 000) | as above | p. 144 | `binding_constant` |
| Klotz binding constant K, 2-pentylpyridine × 11S glycinin | **228 ± 30** | M^-1 (MW 320 000) | as above | p. 144 | `binding_constant` |
| primary binding sites n, SPI / 7S / 11S | **3 / 7.5 / 9** | mol per mol protein | as above; **no uncertainty printed** | p. 144 | `binding_constant` (dimensionless companion) |
| **per-gram binding constant, 2-pentylpyridine × soy protein isolate** | **3.21e-3** | L/g protein | 25 C, water, 5 g/L, equilibrium dialysis | nK/100 000 (mine) | **`derived_assumption`** — the registry's own K_g form, on **the same 100 000 g/mol basis as the shipped Damodaran soy rows** |
| per-gram binding constant, 7S / 11S | **6.14e-3 / 6.41e-3** | L/g protein | as above | nK/MW (mine) | `derived_assumption` |
| **bound amount at 74 °C, SPI / 7S / 11S** | **0.396 ± 0.005 / 0.469 ± 0.006 / 0.569 ± 0.005** | mol 2-pp per mol protein | **74 C**, water, 5 g/L, 0.03 mM 2-pp, 14-20 h dialysis | Table 1, p. 143 | **`measured_ratio`** — the only aqueous binding measurement above 60 C in this batch |
| bound amount at 25 °C, SPI / 7S / 11S | **0.159 ± 0.003 / 0.235 ± 0.003 / 0.375 ± 0.003** | mol per mol protein | 25 C, otherwise as above | Table 1, p. 143 | `measured_ratio` |
| bound amount at 4 °C, SPI / 7S / 11S | **0.241 ± 0.005 / 0.310 ± 0.003 / 0.417 ± 0.005** | mol per mol protein | 4 C, otherwise as above | Table 1, p. 143 | `measured_ratio` |
| **TEMPERATURE FACTOR, 74 C over 25 C** | **2.49 (SPI) / 2.00 (7S) / 1.52 (11S)** | × | water, 5 g/L, 0.03 mM | Table 1 (mine) | **`within_study_ratio`** — and in the linear regime this IS a ratio of binding constants (section 3 item 1) |
| temperature factor, 4 C over 25 C | **1.52 / 1.32 / 1.11** | × | as above | Table 1 (mine) | `within_study_ratio` — **binding is NON-MONOTONE in temperature; 25 C is the minimum** |
| implied K_g at 74 C | **8.00e-3 (SPI) / 1.23e-2 (7S) / 9.73e-3 (11S)** | L/g protein | 74 C, water | K_g(25 C) × the 74/25 ratio (mine) | **`derived_assumption` — DO NOT SHIP** (both steps compromised; Flags 1 and 6) |
| **pH FACTOR, pH 9 over pH 4.5** | **2.76 / 1.93 / 1.76** | × | 25 C, **1 M phosphate buffer** | Table 1 (mine) | **`within_study_ratio`** — 2-pp is a base; protonation releases it |
| bound amount at pH 4.5 / 7 / 9, SPI | **0.119 / 0.141 / 0.328** | mol per mol protein | 25 C, 1 M phosphate | Table 1, p. 143 | `measured_ratio` |
| bound amount at pH 4.5 / 7 / 9, 7S | **0.208 / 0.284 / 0.401** | mol per mol protein | as above | Table 1, p. 143 | `measured_ratio` |
| bound amount at pH 4.5 / 7 / 9, 11S | **0.297 / 0.366 / 0.522** | mol per mol protein | as above | Table 1, p. 143 | `measured_ratio` |
| **IONIC-STRENGTH FACTOR, 1.711 M over 0.171 M NaCl** | **0.22 / 0.41 / 0.54** | × | 25 C | Table 1 (mine) | **`within_study_ratio`** — the registry has NO ionic-strength term |
| bound amount at 0.171 / 0.855 / 1.711 M NaCl, SPI | **0.217 / 0.105 / 0.048** | mol per mol protein | 25 C | Table 1, p. 143 | `measured_ratio` |
| bound amount at three NaCl levels, 7S / 11S | **0.310 / 0.203 / 0.126** and **0.385 / 0.299 / 0.206** | mol per mol protein | 25 C | Table 1, p. 143 | `measured_ratio` |
| **UV FACTOR, exposed over control** | **2.29 / 1.90 / 1.42** | × | 25 C, dual-wavelength 365 + 254 nm during the whole dialysis | Table 1 (mine) | `within_study_ratio` — **and a mass-balance check rules out synthesis for THIS block only** |
| bound amount, UV control / exposure, SPI | **0.154 ± 0.003 / 0.353 ± 0.005** | mol per mol protein | 25 C | Table 1, p. 143 | `measured_ratio` |
| fraction ordering on K | 11S / 7S = **1.74x**, 11S / SPI = **2.13x** | × | 25 C | p. 144 (mine); **the paper says "2-fold"** | `within_study_ratio` |
| the same on a per-gram basis | 7S and 11S agree to **4 %**; whole SPI is **half** of both | × | 25 C | section 3 item 2 (mine) | `within_study_ratio` — **the ordering depends on the basis** |
| saturation concentration, 11S / 7S / SPI | **1.8 mM / 2.4 mM / none to 3 mM** | mM 2-pp | 25 C | p. 144 | `measured_bound` |
| protein content, 7S / 11S fractions | **97.7 % / 92.1 %** protein; **70.4 % / 90.6 %** of the target protein | % dry basis | — | p. 142 | `level_only` |
| protein content, commercial Supro 500E | **not stated** | — | — | — | **absent** |
| **2-pentylpyridine formation pathway** | ammonia + 2,4-decadienal → Schiff base → ring closure → dihydropyridine → oxidation | — | in SPI processing | p. 142 | **`structural_gate`** — `[C]`, Buttery 1977, isotopically confirmed by Zhou & Boatright 2000 |
| amino acids that raise 2-pp / have no effect | **Arg, Lys, Asn, Gln** / **Asp, Glu, Gly, His** | — | in SPI processing | p. 142 | `structural_gate` — **`[C]`, not measured here** |
| odour description | "penetrating-grassy aroma"; "throat-catching taste in water" | — | GC-O and sensory panel | p. 142 | `level_only` — **`[C]`, Boatright & Crum 1997; no odour threshold is printed** |

### Can these be put on the same basis as the shipped binding constants, i.e. converted to K_g in L/g?

**Yes for the 25 C constants, and better than for any other paper in this batch — the basis match with the
registry's existing soy rows is exact. No for the 74 C measurement, which is not a constant.**

- **The molar mass is printed and it is the registry's own soy convention.** Zhou uses **100 000 g/mol for
  whole SPI**, citing O'Neill & Kinsella 1987. Damodaran stated 100 000 g/mol and `kg_2_heptanone_soy`,
  `kg_2_octanone_soy` and `kg_2_nonanone_soy` all carry `"molar_basis": "stated_by_source"` on that figure.
  **So K_g = nK/100 000 = 3.21e-3 L/g for 2-pentylpyridine is directly comparable to 4.40e-3 L/g for
  2-heptanone with no convention mismatch at all.** (Contrast Guo 2019, which uses 220 000 for the same
  protein — a 2.2x conflict.) The 7S and 11S masses (160 000 and 320 000) are also printed and also cited.
- **The protein loading is printed**: 0.5 % w/w = 5 g/L in the protein compartment. **But it is a nominal
  slurry concentration and the paper also measured protein by Bradford** — and never reconciles the two, nor
  states the protein content of Supro 500E, so the *protein* loading may be below 5 g/L by however much the
  isolate is not protein (typically 5-10 %). And the equal-volume dialysis geometry leaves it ambiguous whether
  v is per mole of protein in its own compartment or averaged over the whole system (Flags 4). **A K_g computed
  as nK/MW does not depend on the loading at all** — the loading only matters if one tries to work from Table 1's
  bound amounts instead.
- **The method is `equilibrium_dialysis`**, with provenance qualifications: SPME headspace quantitation via a
  within-run ratio (which cancels calibration offsets, the Amendment 4 logic), and **no thiol-blocking agent**,
  unlike Damodaran's 10 mM 2-mercaptoethanol. For a pyridine with no carbonyl this is a small concern; it would
  be a large one for an aldehyde.
- **The pH is a real, buffered 7.0 for the pH block** (1 M phosphate) and is **unstated for the temperature,
  NaCl and UV blocks**, which appear to be in water. **The 25 C binding constants come from binding curves whose
  medium is never stated** — the Methods describe the buffer only under "For determining pH effects". So a
  `ph_of_measurement` for the 107 M^-1 row would have to be recorded as **unknown**, not 7.0 (Flags 2).
- **The 74 C measurement is a bound amount, not a constant, and cannot be converted.** It has one ligand
  concentration and no isotherm. Section 3 item 5 shows what happens if one chains it through the 25 C
  constant, and why the result should not be shipped.

**On the B26 temperature limit, stated exactly.** This paper **raises the highest aqueous temperature at which
a protein-flavour binding measurement exists in this corpus from 40 C to 74 C** — and it does so on soy, not
pea. What it supplies at 74 C is a **relative** quantity (bound amount against a 25 C control on the same
protein, ratio 1.5-2.5x) and not an absolute constant. **It therefore licenses one specific thing: the
statement that aqueous protein binding does NOT decay with temperature the way a dry-phase or a
physical-adsorption picture would predict, and instead increases by roughly 1.5-2.5x from 25 to 74 C on soy.**
It does not license a pea constant at 90 or 140 C, it does not license a soy constant at 140 C, and it does not
license an absolute soy constant at 74 C. **It does mean that any conservative treatment which assumed binding
falls away at process temperature is measurably wrong in direction on soy.**

**Nothing here goes to `matrix_sites.py`.** No rate constant, no time series (a single 14-20 h equilibration),
no activation energy, and — worth saying explicitly because it is the confusion this repository has flagged
before — **no ΔH, ΔS or ΔG of any kind appears in this paper**, so there is no enthalpy here to misread as an
activation energy. There is also no covalent chemistry: 2-pentylpyridine has no carbonyl, is not a Michael
acceptor, and the paper states its binding "is reversible" (p. 143, as the reason the solvent-extraction route
failed).

## 5. Flags

1. **THE 74 C ARM HOLDS A SOY SLURRY AT 74 C FOR 14-20 HOURS WITH 2-PENTYLPYRIDINE'S OWN PRECURSORS PRESENT,
   AND NO SYNTHESIS CONTROL WAS RUN FOR IT.** This is the most serious problem in the paper and it sits under
   the one result the repository most wants. The same authors established (p. 142, Zhou & Boatright 2000, with
   13C and 15N labelling) that 2-pentylpyridine **forms** in soy protein isolate from ammonia and
   2,4-decadienal, that pro-oxidants and UV increase it, and that arginine, lysine, asparagine and glutamine
   promote it. A soy protein slurry held at **74 C for 14-20 hours** is a set of conditions under which that
   chemistry runs. **The paper ran exactly the right check — a headspace mass balance — and ran it only for the
   UV block**: *"Because the overall 2-pp contents in the headspaces for samples exposed to UV light in this
   equilibrium study were similar to the control, it is unlikely that the increased binding observed in this
   study was the result of 2-pp synthesis"* (p. 144). **No equivalent statement is made for the 74 C block.**
   The isotope-dilution assay does not rescue this: newly synthesised 2-pp is unlabelled and is counted exactly
   like the dosed 2-pp. **Newly formed 2-pp appearing preferentially on the protein side of the membrane would
   register as increased binding.** Until that check exists, the 2.49x rise at 74 C on whole SPI is
   **binding OR synthesis, and the paper cannot separate them.** Note that the fraction with the smallest rise
   (11S, 1.52x) is also the purest and least processed preparation, which is the direction a
   synthesis artefact would predict.
2. **The four treatment blocks do not share a medium, and their 25 C controls disagree by 1.54x.** The pH block
   is in **1 M phosphate buffer**; the temperature, NaCl and UV blocks appear to be in water (no buffer is
   mentioned). The four nominally comparable whole-SPI reference values are **0.141 (pH 7 in 1 M phosphate),
   0.159 (25 C), 0.217 (0.171 M NaCl) and 0.154 (UV control)** — a **1.54x spread (mine)** among conditions that
   should all be "SPI at 25 C". **The pH block is also self-confounded**: 1 M phosphate is a very high ionic
   strength, and this paper's own NaCl block shows that raising ionic strength suppresses binding by up to
   4.5x. So the pH series was measured under conditions the paper itself shows are suppressive, and the pH
   effect and the ionic-strength effect are not separated. **And the medium of the 0.3-3.0 mM binding curves —
   the source of the only binding constants in the paper — is never stated at all.**
3. **The abstract distinguishes "more bound" from "greater affinity" and the paper prints only one quantity.**
   *"More 2-pp was also bound at high temperature (74 °C) than at 25 °C, but greater binding affinity of 2-pp
   was observed at 4 °C than at 25 °C."* Table 1 has a single column type — moles bound per mole of protein —
   and **its 4 C value (0.241) is higher than its 25 C value (0.159) in exactly the same way its 74 C value
   (0.396) is.** There is no separate affinity measurement at 4 C, no binding curve at 4 C or 74 C, and no K at
   any temperature other than 25 C. **The "affinity" language has no printed number behind it and should not be
   quoted as if it did.**
4. **The volume bookkeeping is ambiguous and it moves every derived number.** Equal volumes of protein slurry
   (outside) and 0.03 mM ligand (inside) were used, and the ligand equilibrates across the membrane. **After
   equilibration the ligand averaged over the whole system is 0.015 mM and the protein averaged over the whole
   system is 2.5 g/L.** Table 1's footnote says "the initial concentration of 2-pp in protein slurry was
   0.03 mM", which is ambiguous between the dosing concentration and the post-equilibration concentration in the
   protein compartment, and the paper never says which volume v is normalised on. **A factor-of-two ambiguity
   sits under every bound amount in Table 1.**
5. **The bound amount is inferred from a ratio of two headspace measurements across a membrane, and the two
   compartments are not alike.** B = T(I − O)/I with I and O the headspace concentrations inside and outside.
   The inside compartment is aqueous ligand; the outside contains **5 g/L of a surface-active protein**. Soy
   protein foams, changes surface tension and can itself alter the air/liquid partition of a volatile
   independently of binding. **Any such effect is indistinguishable from binding in this equation.** The paper
   notes elsewhere that dynamic headspace "could not prevent foaming of the protein solutions" (p. 143), so it
   is aware the protein is surface-active. Separate calibration curves were made for the inside and the outside
   — which mitigates a calibration difference but not a partition difference.
6. **The printed binding constants do not reproduce the paper's own Table 1, and are 6-17x too small to do
   so.** Detailed in section 3 item 3. The Klotz fits are on the **0.3-3.0 mM** branch while Table 1 is at
   **0.03 mM** — a factor of ten below the lowest fitted point, i.e. an extrapolation, not an interpolation.
   The most likely explanation is a high-affinity class at low concentration that the fit does not resolve —
   exactly what Guo 2019 found on soy for esters, with a primary class 10-70x above the secondary one. **A
   second internal inconsistency points the same way: Figure 1's y-axis runs to 30 mole/mole protein and all
   three curves rise well above the printed n of 3, 7.5 and 9** (section 3 item 4). **The paper reconciles
   neither.** Whichever branch is meant, **the two candidate per-gram constants for 2-pentylpyridine on whole
   soy differ by 16x: 3.21e-3 L/g from the fitted K, or up to 5.30e-2 L/g implied by Table 1 (mine).**
7. **The 7S "fraction" is only 70.4 % beta-conglycinin.** The paper is admirably explicit — "The purity was such
   that 90.6 % of the 11S fraction was glycinin, and 70.4 % of the 7S fraction was beta-conglycinin" — but it
   means the 7S constant (131 ± 16 M^-1) is a weighted average over a preparation nearly a third of which is
   something else, plausibly including glycinin, which binds 1.74x more strongly. **The true beta-conglycinin
   constant is therefore lower than 131**, and the 1.74x glycinin/beta-conglycinin contrast is an underestimate.
   The 11S value is on firmer ground at 90.6 % purity.
8. **The commercial isolate has no stated protein content.** Supro 500E's protein content, ash, moisture and
   residual lipid are absent, while the two purified fractions have theirs printed. So the one preparation most
   likely to be compared with the registry's other soy rows — a commercial isolate, like Damodaran's and like
   Aspelund's Edi-Pro A — is the one whose composition is unknown. Protein was measured by **Bradford against
   BSA**, under the stated assumption "that all protein constituents have equal binding affinity for Coomassie
   brilliant blue" — an assumption known to fail across soy fractions, and one that propagates directly into v.
9. **Methanol is present in every sample and is never accounted for.** The ligand was dosed as a "2-pp methanol
   solution". The final methanol fraction is never stated, no methanol-only control is run, and methanol is not
   mentioned again. (The same defect appears in Bi 2022 with methanol at ~1.25 % v/v and in Guo 2019 with
   propylene glycol.)
10. **Every molar mass is cited, and the whole-SPI one is called "hypothetical" by the authors.** p. 144: *"the
    data of binding properties calculated on the basis of the hypothetical molecular weight of 100,000 for whole
    SPI"*. It is the right convention for matching the registry's Damodaran rows, and it is still a convention.
    A K in M^-1 and a K_g in L/g both scale inversely with it. The 11S figure of 320 000 is the hexamer and the
    7S figure of 160 000 the trimer, so the two fractions are on structurally sensible bases while the isolate
    is on a round number.
11. **The 4 C result is unexplained and is explained by a citation that says something different.** The paper
    attributes the 4 C rise to Damodaran & Kinsella's observation of increased hydrophobicity at 5 C versus
    20 C, and suggests "the lower temperature allowed protein unfolding". **Cold-induced unfolding of soy
    globulins between 25 and 4 C is not a standard result**, and no structural measurement was made here at any
    temperature. The result stands as measured; the mechanism offered is borrowed and thin. It also means
    **binding is non-monotone in temperature with a minimum at 25 C**, which no simple thermodynamic form
    reproduces and which any interpolation between the three points would miss.
12. **The n values carry no uncertainty and the "2-fold" claim does not match the printed constants.** n = 3,
    7.5 and 9 are printed bare while K carries ± 10, 16 and 30. And the text's *"the glycinin fraction showed a
    2-fold greater affinity for 2-pp than the beta-conglycinin fraction"* is **1.74x** on the printed numbers
    (228/131, mine) — within the combined error bars, but the rounded claim overstates it.
13. **The letters in Table 1 make every difference significant, which is worth a second look.** All twelve
    within-block comparisons on all three proteins carry distinct letters at P < 0.05, on n = 3, with standard
    deviations of 1-2 % of the mean. Standard deviations that tight on a triplicate dialysis-plus-SPME
    measurement are unusual; the paper says 2-pp "was quantified in triplicate by GCD system" (p. 143), which
    reads as **three GC injections**, not three independent dialyses. **If the replication is analytical rather
    than biological, the error bars understate the true variance and the universal significance is an
    artefact.** The paper does not say which it is.
14. **The 74 C dialysis has no membrane-integrity control at that temperature.** "No protein was found inside
    the membrane tube after dialysis" is stated once, without a temperature. A 3 500 MWCO regenerated-cellulose
    membrane held 14-20 h at 74 C is being used well outside its comfortable range, and soy protein at 74 C
    partially denatures and can aggregate or fragment. **If any peptide crossed at 74 C, the outside/inside
    headspace ratio would shift in the direction reported.**
15. **No DOI is printed** (section 0). Cite by volume/issue/page and by the manuscript number MS 20001537.
16. **What this paper does NOT contain**: any binding constant at any temperature but 25 C; any thermodynamic
    quantity; any rate constant; any homologous series (so **no chain-length check**); any aldehyde, ketone,
    ester or 2-alkenal; any covalent chemistry; any odour threshold; any protein content for the commercial
    isolate; any thiol, disulfide or amine assay; any structural measurement (no fluorescence, no CD, no
    hydrophobicity, no particle size) at any temperature; any stated medium for the binding curves; any
    supplementary material.
17. **What to request** (the corresponding author's 2001 address is a University of Kentucky one; recorded for
    completeness): (i) **a 2-pp mass balance for the 74 C arm**, which is the one measurement that would make
    the above-60 C result usable; (ii) binding curves at 4 C and 74 C, which would turn the bound amounts into
    constants; (iii) the medium and pH of the 25 C binding curves; (iv) the reconciliation between the fitted K
    and Table 1's bound amounts at 0.03 mM; (v) whether the triplicate is analytical or independent; (vi) the
    protein content of the Supro 500E lot; (vii) the final methanol fraction.
