# Crowther 1980 — EXTRACTION (Edi-Pro A soy protein isolate autoclaved 20 min at 100 or 121 C at 29 % or 40 % moisture, or extruded at 71 C, then dried, ground and packed dry into a GC column; adsorption coefficients K in mL/m^2 at 70 C and heats of adsorption from 60/70/80 C, for three 2-alkanones, three 1-alkanols and hexanal)

### THE PREHEAT EXPERIMENT `matrix_sites.py` SAYS IT DOES NOT MODEL. That module charges its binding sites once at the start of the cook and does not change them with heating. This paper heats soy protein isolate at **100 and 121 C** — real process temperatures — and then measures what its flavour binding has become, against an unheated control from the same product. The answer is a clean, quantified split: **the adsorption coefficient falls by 47-49 % for the ketones, 35-40 % for the alcohols and 43 % for hexanal (Table 6, p. 109), while the heat of adsorption does not change at all (Table 5, p. 108: every ΔH p-value between 0.076 and 0.968).** The authors' own reading is the one the repository would need: *"Since ΔH values were unaffected by processing, the observed changes in K are postulated to be the result of changes in the NUMBER OF BINDING SITES"* (p. 112). **Heating removes sites; it does not weaken them.** And a second arm says the sign is not automatic — **shear at 71 C RAISES alcohol adsorption by 14-46 %** (p. 111). Everything here is measured on **dry** protein by inverse gas chromatography, so **no row is shippable into `REVERSIBLE_BINDING`**; what is shippable in spirit is the *within-study ratio*, and it is the only measurement of a preheat effect on flavour binding in the corpus that has an unheated control on the same lot.

**Source on disk:** `data/articles/crowther1980.pdf` (17 pp., J. Food Process Eng. **4** (1980) 99-115).
Read from the `pdftotext -layout` text layer, with **Tables 2 (p. 105), 3 (p. 106), 4 (p. 107), 5 (p. 108),
6 (p. 109) and 7 (p. 111) all verified cell-by-cell against the rendered page images** — necessary, because
Tables 3-6 are printed **rotated 90 degrees** and the text layer scrambled or mangled them badly (Table 3's
fifth column header came out as "Hexanal" when the page reads **Hexanol**; the untreated hexanone K came out as
3.45 when the page reads **3.48**; the 121 C / 29 % heptanone standard error came out as "$6.63" when the page
reads **±0.03**; the 121 C / 40 % hexanal K came out as 1.78 when the page reads **1.79**; Table 4's first
pentanone cell came out as "+.59" when the page reads **−6.59**; the sheared pentanone ΔH came out as −9.12 when
the page reads **−9.11**; Table 7's lactose hexanol came out as 17.67 when the page reads **17.57**). The
re-typed tables below are from the page images. **Figures 1-4 are images**: Fig. 1 shows nine representative
chromatographic peaks with no printed numbers, Figs. 2 and 3 are the ln(V_m/T) vs 1/T regression lines for the
alcohols and the carbonyls on one 121 C / 29 % packing and carry only **correlation coefficients** inside the
panels (heptanone 0.9930, hexanone 0.9973, hexanal 0.9978, pentanone 0.9378 — legible on Fig. 3; the Fig. 2
alcohol legend is not legible), and Fig. 4 is four scanning electron micrographs. **Adsorption coefficients at
60 C and 80 C are stated to exist and are explicitly NOT SHOWN** (p. 103): "The complete data are available
(Crowther 1979)" — an unpublished thesis, not on disk. Repo status before this dossier: Crowther 1980 has **no
extraction dossier** and is **not cited** in `src/kinetic_core/parameters_matrix.py`,
`src/kinetic_core/matrix_sites.py` or `data/species/protein_matrices.yml`.

## 0. Identity

| field | value |
|---|---|
| Title | "Effects of Processing on Adsorption of Off-Flavors onto Soy Protein" |
| Authors | **Alan Crowther** (Celanese Chemical Company, Bishop, Texas — the work was done at Iowa State on a Celanese Corporation fellowship), **Lester A. Wilson** (Department of Food Technology and Agricultural Experiment Station, Iowa State University, Ames, Iowa), **Charles E. Glatz** (Department of Chemical Engineering and Engineering Research Institute, Iowa State University) |
| Venue | **Journal of Food Process Engineering 4 (1980) 99-115** — exactly as printed in the footer of p. 99. "© Copyright 1981 by Food & Nutrition Press, Inc., Westport, Connecticut" |
| DOI | **NO DOI IS PRINTED IN THE PDF.** Nothing DOI-shaped appears on any of the 17 pages. |
| Dates | **Received for Publication October 24, 1980; Accepted for Publication January 15, 1981.** Portions presented at the IFT Annual Meeting, June 10, 1980, New Orleans |
| **Year ambiguity, flagged** | The journal line says volume 4 **(1980)**, the copyright line says **1981**, and Aspelund & Wilson 1983 — the same laboratory — cite this paper as **"Crowther, A.; Wilson, L. A.; Glatz, C. E. *J. Food Process Eng.* 1981, *4*, 99."** The filename and this dossier use 1980 because that is what the journal footer prints. **A citation search on "Crowther 1981" and on "Crowther 1980" retrieves the same paper.** |
| Funding | Iowa State University Research Foundation; Iowa Agriculture and Home Economics Experiment Station (Paper No. J-9920, Project No. 2164, contributing to North Central Regional Research Project NC-136); Engineering Research Institute of Iowa State University; A.C. on a Celanese Corporation fellowship |
| Protein | **Edi-Pro A**, isoelectric soy protein isolate, Ralston-Purina Company, St. Louis, Missouri — **"reported by the company to contain 93.5 % protein, 5.5 % moisture, and 0.2 % fat"** (p. 100). This is the **same commercial product Aspelund 1983 used**, and unlike Aspelund this paper prints its composition. |
| The 8 ligands | n-hexane and 2-pentanone (Matheson, Coleman & Bell); 2-hexanone, 2-heptanone, n-hexanal, 1-pentanol (Aldrich); 1-hexanol (Eastman Kodak); 1-heptanol (J. T. Baker). **n-hexane "was not significantly retained" and drops out entirely** (p. 103), leaving 7 |
| Naming | "K" = adsorption coefficient, **mL/m^2**, a Henry's-law constant C_s/C_g at 70 C; "ΔH" = heat of adsorption, kcal/g-mol, from the slope of ln(V_m/T) vs 1/T; "V_m" = normalised retention volume, mL/g; "V_R" = corrected retention volume, mL; "PDI" = protein dispersibility index; "S" = specific surface area, m^2/g |
| Companions on disk | `aspelund1983_extraction.md` (**the direct predecessor: same laboratory, same Edi-Pro A, same method, untreated protein at 80/90/100 C — and the two papers disagree, see Flags 4**), `guo2019_extraction.md` and `Xu2022_extraction.md` (the modern aqueous preheat studies on soy and pea), `damodaran1981_extraction.md`, `bi2022_extraction.md` |

## 1. Why it matters

**It is the preheat experiment, and the repository has an explicit hole where it goes.**
`src/kinetic_core/matrix_sites.py` charges its reactive pools — free thiol, disulfide, amine, in mmol per gram
of protein — **once, from `data/species/protein_matrices.yml`, at the start of the thermal programme**, and then
integrates aldehyde and HMF binding over the cook against those fixed pools. Whether heating the protein itself
changes the number of available sites is not modelled. This paper is a direct, controlled test of exactly that
question, on soy, with an unheated control from the same commercial lot:

- **Autoclaving at 100-121 C cuts the adsorption coefficient by 35-49 %**, and the paper's own summary of which
  by how much is printed on p. 111: *"autoclaving greatly reduces adsorption (47-49 % for ketones, 35-40 % for
  alcohols, and 43 % for hexanal)"*. Recomputed from Table 3 (mine): pentanone **−47.1 %**, hexanone **−48.9 %**,
  heptanone **−49.0 %**, pentanol **−34.5 %**, hexanol **−35.0 %**, heptanol **−39.9 %**, hexanal **−43.0 %**.
- **The heat of adsorption does not move.** Table 5's ANOVA gives, for ΔH, p-values of 0.076-0.974 across every
  compound and every factor — nothing significant anywhere — and the text says the t-tests were not even run on
  ΔH "since the data clearly indicated that there was no significant change in this parameter" (p. 103).
- **So the mechanism is site COUNT, not site STRENGTH.** p. 112: *"Since ΔH values were unaffected by
  processing, the observed changes in K are postulated to be the result of changes in the number of binding
  sites."* And the proposed physical picture, p. 113: *"a chain unfolding of the globular proteins, which
  exposes the more nonpolar regions of the protein, decreasing both solubility and the availability of polar
  hydrogen bonding sites."* **This maps directly onto the `matrix_sites.py` data model**, whose site densities
  are counts per gram and whose rate brackets are the strength term. It says: on heating, scale the counts, not
  the rate.
- **And the corroborating denaturation gradient is in the same table.** Table 2 (p. 105) prints protein
  dispersibility index for the same samples: neutralised PDI **88.4 (untreated) → 48.7 (sheared) → 7.8
  (autoclaved)**, with disc gel electrophoresis going from "several strong bands" to "bands less distinct ...
  and some missing" to **"No bands"**. The protein that lost 90 % of its dispersibility lost 35-49 % of its
  adsorption. **This is the corpus's only paired denaturation-and-binding measurement.**

**But the sign is not automatic, and that is the second reason to read it.** The extruded (sheared, 71 C
maximum barrel temperature, 33 % moisture) sample went the *other* way for the alcohols: **pentanol +46 %,
hexanol +37 %, heptanol +14 %** against untreated (p. 111, and recomputed from Table 3, mine), while its
ketones were flat to slightly down (pentanone +4.5 %, hexanone **−6.9 %**, heptanone −3.4 %) and hexanal +9.0 %.
A model that assumed "processing reduces binding" as a rule would get the sheared alcohols wrong by up to 1.5x
in the wrong direction. The registry has been here before: two of Leksrisompong's three caseinate rows are
**negative** — the matrix makes the odourant *more* volatile — and `parameters_matrix.py` carries them precisely
so the layer "can emit a shift below 1 at all". **Crowther's sheared alcohols are the preheat analogue of that
lesson.**

**The B26 temperature limit: what this paper does and does not do to it.** Wave B26's record says "37 C is an
in-mouth temperature, not a process one: nothing here licenses a pea binding constant at 90 or 140 C." Two
temperatures matter here and they are different things:

- **The TREATMENT temperature is 100 and 121 C** — genuinely inside the process band, held 20 min in a preheated
  autoclave at controlled moisture. This is the first source in the corpus that puts a protein through a real
  thermal process and *then* measures its flavour binding.
- **The MEASUREMENT temperature is 60, 70 and 80 C**, with all reported K values at **70 C**. That too is above
  the 60 C line the B26 note draws, and it is on soy.
- **But the measurement is dry gas-solid inverse GC**, exactly as in Aspelund 1983, so it still does not supply
  an *aqueous* constant above 60 C. **The limit is not lifted by this paper either.** What it supplies is a
  ratio — how much a heat treatment changes binding — and a ratio is the object `parameters_matrix.py` was
  built to carry (Amendment 4).

**The chain-length slope, checked a third time.** From the untreated K values at 70 C (Table 3, mine): the
2-alkanone series gives **2.23x and 1.93x per CH2** (geometric mean **2.08x**) and the 1-alkanol series **1.92x
and 2.40x** (geometric mean **2.15x**). Against the shipped `CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81` that is **24-26 %
low**, in the same direction as Aspelund's 2.23-2.27x at 90 C. **And the slope barely moves with treatment**:
autoclaved ketones **2.04x**, sheared ketones **2.00x** (mine) against untreated 2.08x. That is a genuinely
useful null — *whatever heating does to soy binding, it does not change the chain-length dependence*, so a
chain-length transfer rule survives a preheat where an absolute constant does not.

**The functional-group ordering is again the dry-phase one, and the paper says so in the sharpest terms in the
corpus.** p. 112: *"The comparative extents of binding (as indicated by K) in this dry system are predominantly
**the reverse of the order found in aqueous solutions of soy proteins** (Beyeler and Solms 1974; Gremli 1974).
Here, alcohols clearly had the highest K's while the carbonyls were quite close."* At six carbons, untreated,
70 C: hexanol **7.83x** hexanal and **6.24x** hexanone (mine). In aqueous soy the ordering Beyeler & Solms
report, quoted on p. 100, is **aldehyde > ketone > alcohol** — the exact inverse. **This is direct measured
support for `parameters_matrix.py`'s refusal of any log-P-shaped matrix term (k4b hold-out guard #4), and a
hard prohibition on transferring any functional-group contrast from this paper to an aqueous model.**

**One more thing the paper contributes: an explicit reversible/irreversible separation.** p. 112: *"In this
study, the only molecules detected were those that traversed the column and thus were only physically
(reversibly) adsorbed"* — the method **cannot see** an irreversibly bound molecule, because such a molecule
never elutes. That is a genuinely useful property. Contrast Bi 2022, whose headspace-depletion binding percentage
counts irreversible capture as binding and needed quarantining for (E)-2-octenal on exactly that ground. **An
inverse-GC K is, by construction, a reversible-only quantity.** The paper reasons about the difference explicitly
against Arai 1970 (which found *more* hexanal than hexanol retained on denatured soy — the opposite of the K
ordering here) and concludes Arai measured chemical, irreversible bonding while this study measures physical.
`kg_hexanal_soy_denatured` (1.47e-3 L/g) in the registry comes from **Arai 1970 via Damodaran** — the very study
Crowther argues was measuring a different phenomenon (Flags 7).

**What this paper does NOT give the repository**: any aqueous measurement; any binding constant in M^-1 or L/g;
any pH; any protein concentration in g/L; any rate constant or activation energy; any covalent adduct
measurement; any water activity (the authors say so themselves and call it the variable they should have
measured); the 60 C and 80 C adsorption coefficients (stated to exist, not shown); any 2-alkenal or
alpha,beta-unsaturated carbonyl (so **nothing touches `ALPHA_BETA_UNSATURATION`**); any pyrazine, pyridine,
furan or sulfur compound; any error bar on the untreated or sheared ΔH values.

## 2. Methods as they matter to a model

- **The protein and its composition.** Edi-Pro A, isoelectric soy protein isolate (Ralston-Purina), **93.5 %
  protein, 5.5 % moisture, 0.2 % fat as reported by the manufacturer** (p. 100). No thiol, disulfide or amine
  assay; no molar mass; no SDS-PAGE beyond the qualitative disc gel of Table 2.
- **The treatments — this is the experiment.** A **2 × 2 factorial**: temperature (**100 C, 121 C**) ×
  moisture content (**29 %, 40 %**), plus a separate shear treatment and an untreated control. **All treatments
  repeated twice to provide replicate samples** (p. 101).
  - **Moisture** was set by adding water to the protein to a homogeneous mixture, then **drying to the target
    weight in an air-circulation drier at 55 C with occasional stirring**.
  - **Heat**: the protein was **spread thinly on trays and heated 20 min in a PREHEATED AUTOCLAVE** at 100 or
    121 C. So the thermal load is 20 min at a plateau, moist, in saturated steam conditions — a retort process,
    not a dry bake.
  - **Shear**: isolate at **33 % moisture** fed to a **Wenger X-5 laboratory extruder at 700 rpm without a die**
    (hence little back pressure). **Maximum barrel temperature reached was 71 C.** So the shear arm is a low-heat,
    high-mechanical-energy treatment and is not comparable to the autoclave arm on temperature.
- **Return to a common state before measurement.** Every sample — treated and untreated — was **dried back to
  the initial Edi-Pro A moisture content** in an air-circulation drier, **ground in a cyclone grinder**, and
  **sized to 230/325 mesh with Tyler sieves** (p. 101). The sizing is important: it puts every packing in the
  same particle-size window, which is why the surface areas come out within 0.17-0.28 m^2/g across treatments
  despite very different processing.
- **Column and chromatography.** **3-ft glass columns, 2 mm ID**, cleaned with 40 % NaOH, rinsed with distilled
  water then acetone, dried at 100 C. **Approximately 1.4 g of the treated samples and 0.5 g of unprocessed
  protein** could be packed while keeping reasonable flow at **50 psig head pressure** (p. 101) — note that the
  untreated control was packed at **less than half the mass** of the treated samples (Flags 5). Varian 3700 GC
  with FID and a Varian CDS 111 integrator. Table 1 conditions (p. 102): **injector 150 C, detector 150 C, column
  oven 60 / 70 / 80 C, nitrogen 23 or 28 mL/min, hydrogen 30 mL/min, air 300 mL/min.** Flow rates measured with
  a bubble flowmeter at the detector outlet.
- **Dosing.** "The compounds injected were **drawn from equilibrated reagent headspaces** with gas-tight
  syringes. Where more dilute samples were required, the syringes were flushed several times before injection"
  (p. 101). **No injected mass is quantified anywhere in this paper** — unlike Aspelund 1983, which prints
  10^-9 to 10^-6 g per injection. The dilution is adjusted per compound by an undocumented number of syringe
  flushes, to keep peaks from tailing (Flags 6). The paper asserts the linear (Henry's-law) region is being
  worked in — "Assuming that we are working in the linear portion of the adsorption isotherm (very low gas
  concentrations)" (p. 102) — **as an assumption, not a demonstration**.
- **Measurement family: INVERSE GAS CHROMATOGRAPHY (gas-solid), the same sixth family as Aspelund 1983.** Not
  `headspace_depletion`, not `equilibrium_dialysis`, not `gel_filtration`, not `static_headspace_partition`,
  not `sensory_BET`. **By construction it detects only reversibly adsorbed molecules** (p. 112).
- **The chain of definitions, printed as Eqs. 1-7 (pp. 102-103).**
  - Eq. 1: uncorrected retention volume **V_R' = w(t_R − t_D)**, with w the gas flow rate, t_R the sample
    retention time and t_D the **methane** retention time.
  - Eq. 2: corrected for pressure drop and for the difference between column and flowmeter conditions, using
    **j = (3/2)[(P_i/P_o)^2 − 1] / [(P_i/P_o)^3 − 1]**, the **compressibility correction of James & Martin
    (1952)**, and the column/flowmeter temperature ratio.
  - Eq. 3: normalised on packing mass, **V_m = V_R / m**, in mL/g.
  - Eq. 4: **the heat of adsorption from the slope of ln(V_m / T) against 1/T** (Kiselev & Yashin 1969). Note
    the ordinate is **ln(V_m/T)**, not ln(V_m) and not ln(t_R) — this differs from Aspelund 1983, which
    regressed ln(t_cor) against 1/T with no T in the ordinate (Flags 4).
  - **Eq. 5 is the definition that matters: K = V_m / S**, "an adsorption coefficient independent of packing
    mass or surface area", and *"By definition, K is a Henry's type constant relating the adsorbed surface
    concentration, C_s, to the gas phase concentration in equilibrium with it, C_g: K = C_s/C_g."* **Units:
    mL/m^2** (Nomenclature, p. 114). **The same units and the same construction as Aspelund's V_S**, which is
    what makes the two papers directly comparable — and they do not agree (Flags 4).
  - Eqs. 6-7: **ΔG = −RT ln K** and the entropy from it, "Both parameters will be dependent on the units of K.
    **We consider K and ΔH the more fundamental parameters.**" — so unlike Aspelund, this paper declines to
    report ΔG and ΔS at all, on the correct ground that they are unit-dependent.
- **Surface area.** Measured **for each protein sample** by Micromeritics, **single-point BET with ARGON as
  adsorbate, ±3 % accuracy**; one sample (100 C, 40 % H2O) was done by **multiple-point BET with KRYPTON, ±1 %**
  (Table 2 footnotes, p. 105).
- **Denaturation characterisation.** Scanning electron micrographs (60/40 gold/palladium coating, Fig. 4);
  **disc gel electrophoresis**; and **protein dispersibility index (PDI**, Smith & Circle 1972), each measured
  **both without pH adjustment and after neutralising with 40 % NaOH** (p. 102) — which is why Table 2 has two
  PDI columns, "Isoelectric" and "Neutralized".
- **Statistics.** ANOVA on the 2 × 2 factorial via SAS; because the moisture-temperature **interaction**
  dominated and obscured the main effects (Table 5), the main effects were re-assessed by **Student's t-test**
  (Table 6). **Note the significance coding in Table 6 is inverted from the usual convention** and must be read
  carefully (Flags 3).
- **Replication, stated per row.** Table 3 footnotes: **1** = average of independently autoclaved replicates;
  **2** = average of four entries, two independent replicates each run in two columns; **3** = average of two
  columns packed with material from the same batch. Table 4 footnotes: **1** = in three cases heptanol was not
  detected at 60 C and **only the 70 and 80 C data were used**; **2** = average of two independent replicates;
  **3** = average of four entries, two pairs of columns packed with two independent replicates; **4** = average
  of two columns with the same packing. **The footnote numbering differs between Tables 3 and 4 — the same
  superscript does not mean the same thing in the two tables.**

## 3. Tables re-typed

Every cell below was read off the rendered page image (Tables 3-6 required rotating the page 90 degrees).
Evidence marks: `[M]` measured in this study, `[C]` cited from elsewhere, `[F]` fitted/derived by the authors.

### Table 1 (p. 102). "Experimental chromatographic conditions for the Varian 3740 gas chromatograph"

| condition | value |
|---|---|
| Injector temperature, °C | 150 `[M]` |
| Detector temperature, °C | 150 `[M]` |
| Column oven temperature, °C | 60; 70; 80 `[M]` |
| Gas flow rates, ml/min — Nitrogen | 23, 28 `[M]` |
| Gas flow rates, ml/min — Hydrogen | 30 `[M]` |
| Gas flow rates, ml/min — Air | 300 `[M]` |

The table caption says "Varian 3740" while the Methods text (p. 101) says the analyses were performed on a
**Varian 3700**. One of the two is a typographical slip and the paper does not resolve it.

### Table 2 (p. 105). "Physicochemical properties of the protein samples used as column packings. Temperature and moisture refer to the condition of the sample when autoclaved, not while in the column"

| Treatment | Specific Surface Area^1 (m^2/g) | PDI^2 Isoelectric | PDI^2 Neutralized | Electrophoresis^3 |
|---|---|---:|---:|---|
| **Autoclaved** | | **1.9** `[M]` | **7.8** `[M]` | **No bands** `[M]` |
| 121 °C, 29 % H2O | 0.17, 0.18 `[M]` | | | |
| 121 °C, 40 % H2O | 0.20, 0.24 `[M]` | | | |
| 100 °C, 29 % H2O | 0.23, 0.19 `[M]` | | | |
| 100 °C, 40 % H2O | 0.20, 0.28^4 `[M]` | | | |
| **Untreated** | **0.20** `[M]` | **1.3** `[M]` | **88.4** `[M]` | Several strong bands `[M]` |
| **Sheared** | **0.28** `[M]` | **1.6** `[M]` | **48.7** `[M]` | Bands less distinct than Edi-Pro A and some missing `[M]` |

Footnotes exactly as printed: *^1 Single-point BET ± 3 % accuracy according to Micromeritics Instrument
Corporation. The two values for each autoclaved condition are for the two independently prepared samples.*
*^2 Average of two determinations.* *^3 Neutralized samples.* *^4 This sample was done with a multiple-point BET
analysis (± 1 % accuracy) using krypton as the adsorbate.*

**Read the layout carefully: the PDI and electrophoresis row labelled "Autoclaved" is a SINGLE value for all
four autoclave conditions pooled**, not four values. The four surface-area pairs below it are per-condition.
So there is **no per-condition PDI** and the paper's claim that "both electrophoresis and PDI results changed
most drastically for the moist-heat treated samples" (p. 112) cannot be resolved to which moist-heat treatment.

### Table 3 (p. 106). "Adsorption coefficients for the adsorption at 70 °C of various compounds onto Edi-Pro A subjected to three types of treatment"

Values are the **adsorption coefficient K in mL/m^2 (± standard error of the mean)**, all at **70 °C**.

| Column Packing Treatments | Pentanone | Hexanone | Heptanone | Pentanol | Hexanol | Heptanol | Hexanal |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Autoclaved** | | | | | | | |
| 121 °C, 29 % H2O^1 | 0.822 ± 0.090 `[M]` | 1.89 ± 0.06 `[M]` | 3.49 ± 0.03 `[M]` | 6.66 ± 0.80 `[M]` | 13.5 ± 1.3 `[M]` | 30.4 ± 3.2 `[M]` | 1.54 ± 0.11 `[M]` |
| 121 °C, 40 % H2O^2 | 0.879 ± 0.062 `[M]` | 1.91 ± 0.22 `[M]` | 3.94 ± 0.64 `[M]` | 9.04 ± 0.71 `[M]` | 17.0 ± 2.1 `[M]` | 38.4 ± 4.2 `[M]` | 1.79 ± 0.27 `[M]` |
| 100 °C, 29 % H2O^1 | 1.046 ± 0.023 `[M]` | 2.13 ± 0.20 `[M]` | 3.76 ± 0.50 `[M]` | 8.39 ± 1.65 `[M]` | 16.3 ± 2.8 `[M]` | 35.3 ± 5.1 `[M]` | 1.90 ± 0.30 `[M]` |
| 100 °C, 40 % H2O^1 | 0.556 ± 0.023 `[M]` | 1.20 ± 0.04 `[M]` | 2.62 ± 0.50 `[M]` | 5.51 ± 2.10 `[M]` | 9.6 ± 2.9 `[M]` | 21.2 ± 7.4 `[M]` | 1.09 ± 0.02 `[M]` |
| **Mean** | 0.825 ± 0.102 `[F]` | 1.78 ± 0.201 `[F]` | 3.43 ± 0.32 `[F]` | 7.40 ± 0.81 `[F]` | 14.1 ± 1.7 `[F]` | 31.3 ± 3.8 `[F]` | 1.58 ± 0.18 `[F]` |
| **Untreated^3** | 1.56 ± 0.06 `[M]` | 3.48 ± 0.01 `[M]` | 6.73 ± 0.53 `[M]` | 11.3 ± 0.15 `[M]` | 21.7 ± 1.2 `[M]` | 52.1 ± 1.5 `[M]` | 2.77 ± 0.45 `[M]` |
| **Sheared^3** | 1.63 ± 0.05 `[M]` | 3.24 ± 0.11 `[M]` | 6.50 ± 0.26 `[M]` | 16.5 ± 0.87 `[M]` | 29.7 ± 2.3 `[M]` | 59.4 ± 3.2 `[M]` | 3.02 ± 0.07 `[M]` |

Footnotes exactly as printed: *^1 Average of independently autoclaved replicates.* *^2 Average of four entries.
Two independent replicates each run in two columns.* *^3 Average of two columns packed with material from the
same batch.*

**Internal check (mine): does the "Mean" row equal the mean of the four autoclaved rows?** pentanone
(0.822+0.879+1.046+0.556)/4 = **0.826** against 0.825 ✓; hexanone **1.783** against 1.78 ✓; heptanone **3.453**
against 3.43 (0.7 % off); pentanol **7.400** against 7.40 ✓; hexanol **14.10** against 14.1 ✓; heptanol
**31.325** against 31.3 ✓; hexanal **1.580** against 1.58 ✓. **The Mean row closes on six of seven; heptanone is
0.7 % off, consistent with rounding of the underlying replicates.**

### Table 4 (p. 107). "Heats of adsorption for various compounds onto Edi-Pro A subjected to three types of treatment. Except as noted for heptanol, the ΔH values are based on data at 60, 70 and 80 °C"

Values are the **heat of adsorption in kcal/mol (± standard error of the mean)**. Printed with a negative sign.
Note the column header in the printed table reads **"Hepanol"** — a typographical slip for heptanol.

| Column Packing Treatments | Pentanone | Hexanone | Heptanone | Pentanol | Hexanol | Heptanol^1 | Hexanal |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Autoclaved** | | | | | | | |
| 121 °C, 29 % H2O^2 | −6.59 ± 1.67 `[M]` | −9.31 ± 0.30 `[M]` | −9.86 ± 0.08 `[M]` | −14.7 ± 0.65 `[M]` | −16.6 ± 0.01 `[M]` | −16.2 ± 1.78 `[M]` | −7.92 ± 1.00 `[M]` |
| 121 °C, 40 % H2O^3 | −6.55 ± 1.15 `[M]` | −9.20 ± 0.79 `[M]` | −10.9 ± 0.14 `[M]` | −15.1 ± 0.16 `[M]` | −16.9 ± 0.17 `[M]` | −17.6 ± 0.99 `[M]` | −9.36 ± 1.41 `[M]` |
| 100 °C, 29 % H2O^2 | −6.26 ± 1.10 `[M]` | −8.13 ± 1.14 `[M]` | −10.9 ± 0.11 `[M]` | −16.0 ± 1.75 `[M]` | −17.1 ± 0.28 `[M]` | −17.3 ± 0.06 `[M]` | −8.74 ± 0.01 `[M]` |
| 100 °C, 40 % H2O^2 | −10.0 ± 0.99 `[M]` | −10.7 ± 3.20 `[M]` | −12.7 ± 1.08 `[M]` | −13.9 ± 1.02 `[M]` | −16.4 ± 0.29 `[M]` | −17.8 ± 0.18 `[M]` | −11.4 ± 0.88 `[M]` |
| **Mean** | −7.35 ± 0.89 `[F]` | −9.34 ± 0.53 `[F]` | −11.2 ± 0.6 `[F]` | −14.9 ± 0.4 `[F]` | −16.8 ± 0.2 `[F]` | −17.2 ± 0.4 `[F]` | −9.36 ± 0.74 `[F]` |
| **Untreated^4** | −10.1 `[M]` | −11.0 `[M]` | −12.7 `[M]` | −13.1 `[M]` | −16.7 `[M]` | −16.8 `[M]` | −11.1 `[M]` |
| **Sheared^4** | −9.11 `[M]` | −10.5 `[M]` | −12.3 `[M]` | −14.2 `[M]` | −15.0 `[M]` | −16.3 `[M]` | −10.8 `[M]` |

Footnotes exactly as printed: *^1 In three cases, heptanol was not detected at 60 °C, and only the 70° and 80 °C
data were used.* *^2 Average of two independent replicates.* *^3 Average of four entries. Two pairs of columns
were packed with two independent replicates.* *^4 Average of two columns with the same packing.*

**The untreated and sheared rows carry NO standard errors** — the only rows in either table without them.

**Internal check (mine): the "Mean" row of the four autoclaved values.** pentanone
(6.59+6.55+6.26+10.0)/4 = **7.35** ✓; hexanone **9.335** vs 9.34 ✓; heptanone **11.09** vs 11.2 (1 % off);
pentanol **14.925** vs 14.9 ✓; hexanol **16.75** vs 16.8 ✓; heptanol **17.225** vs 17.2 ✓; hexanal **9.355** vs
9.36 ✓. **Closes.**

**Note the direction of the ΔH's, which the abstract's claim of "no change" partly conceals.** For the
**ketones and hexanal**, every autoclaved value is *less* negative than the untreated one (pentanone −6.26 to
−10.0 against untreated −10.1; hexanal −7.92 to −11.4 against −11.1). For the **alcohols** the autoclaved
values are mostly *more* negative than untreated (pentanol −13.9 to −16.0 against −13.1). **The ANOVA of
Table 5 found none of this significant** — the standard errors are large (up to ±3.20 on one cell) and the
untreated row has no error bar at all to test against.

### Table 5 (p. 108). "Analysis of variance for the 2 × 2 factorial design for effects of temperature and moisture on adsorption coefficients (at 70 °C) and heats of adsorption (from data at 60, 70, and 80 °C). The values are the probability that the independent variable did not affect the dependent variable"

| Dependent Variable | Independent Variable | Pentanone | Hexanone | Heptanone | Pentanol | Hexanol | Heptanol | Hexanal |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **K** | Moisture | 0.107 `[M]` | 0.125 `[M]` | 0.781 `[M]` | 0.819 `[M]` | 0.880 `[M]` | 0.974 `[M]` | 0.509 `[M]` |
| | Temperature | 0.493 `[M]` | 0.237 `[M]` | 0.280 `[M]` | 0.393 `[M]` | 0.290 `[M]` | 0.220 `[M]` | 0.418 `[M]` |
| | Moisture/Temperature Interaction | **0.025** `[M]` | **0.047** `[M]` | 0.177 `[M]` | 0.133 `[M]` | 0.102 `[M]` | 0.103 `[M]` | 0.082 `[M]` |
| **ΔH** | Moisture | 0.414 `[M]` | 0.588 `[M]` | 0.082 `[M]` | 0.556 `[M]` | 0.586 `[M]` | 0.428 `[M]` | 0.197 `[M]` |
| | Temperature | 0.325 `[M]` | 0.929 `[M]` | 0.076 `[M]` | 0.968 `[M]` | 0.901 `[M]` | 0.713 `[M]` | 0.338 `[M]` |
| | Moisture/Temperature Interaction | 0.218 `[M]` | 0.466 `[M]` | 0.174 `[M]` | 0.285 `[M]` | 0.081 `[M]` | 0.569 `[M]` | 0.368 `[M]` |

**Read the caption literally: these are the probabilities that the variable did NOT have an effect**, so a
small number means an effect. On **K**, the only two values below 0.05 are the **moisture × temperature
interaction for pentanone (0.025) and hexanone (0.047)**; every main effect on K is above 0.10. On **ΔH**,
nothing reaches 0.05 anywhere — the smallest are heptanone's temperature (0.076) and moisture (0.082) and
hexanol's interaction (0.081). **This is the entire statistical basis for the paper's central claim that ΔH is
unchanged by processing, and it is a null result on n = 4 conditions** (Flags 2).

### Table 6 (p. 109). "Comparison of adsorption coefficients for various process variables. The values presented are differences between adsorption coefficinets [sic] reported in Table 3."

Values are **ΔK at 70 °C**, in mL/m^2. All are **differences**, not ratios.

| Comparison | Pentanone | Hexanone | Heptanone | Pentanol | Hexanol | Heptanol | Hexanal |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Autoclaved^1** | | | | | | | |
| 29 % H2O | −0.224\*\*\* `[F]` | −0.24 `[F]` | −0.27 `[F]` | −1.73 `[F]` | −2.8 `[F]` | −4.9 `[F]` | −0.36 `[F]` |
| 40 % H2O | 0.323\*\* `[F]` | 0.71\*\* `[F]` | 1.42\*\*\*\* `[F]` | 3.53\*\*\*\* `[F]` | 7.4\*\*\* `[F]` | 17.2\*\*\* `[F]` | 0.7\*\*\* `[F]` |
| 100 °C | −0.49\*\* `[F]` | −0.93\*\* `[F]` | −1.24\*\*\* `[F]` | −2.88 `[F]` | −6.7\*\*\* `[F]` | −14.1\*\*\* `[F]` | −0.81\*\*\* `[F]` |
| 121 °C | 0.057 `[F]` | 0.02 `[F]` | 0.45 `[F]` | 2.38 `[F]` | 3.5 `[F]` | 8.0 `[F]` | 0.25 `[F]` |
| **Untreated vs. Sheared^2** | 0.07 `[F]` | −0.24\*\*\*\* `[F]` | −0.23 `[F]` | 5.2\*\* `[F]` | 8.0\*\* `[F]` | 7.3\*\*\*\* `[F]` | 0.25 `[F]` |
| **Untreated vs. Autoclaved^3** | −0.735\*\* `[F]` | −1.7\* `[F]` | −3.3\* `[F]` | −3.9\*\*\* `[F]` | −7.6\*\*\* `[F]` | −20.8\*\* `[F]` | −1.19\*\*\* `[F]` |

Footnotes exactly as printed: *\*P (probability of no significant difference) < 0.01; \*\* P < 0.05;
\*\*\* P < 0.1; \*\*\*\* P < 0.2.* *^1 (K_121 − K_100) at constant % H2O or (K_40% − K_29%) at constant T;
df = 4 (using an estimate of the standard deviation by pooling the estimates at the four conditions).*
*^2 (K_sheared − K_untreated); df = 2.* *^3 (K̄ − K_untreated) where K̄ is the mean of the four autoclaved
samples; df = 4.*

**Two things must be read exactly or this table will be misused.**
(a) **The star convention is INVERTED from the usual one**: one star is the *strongest* result (P < 0.01) and
four stars the *weakest* (P < 0.2). The two rows with single stars — hexanone and heptanone under
"Untreated vs. Autoclaved" — are the most significant differences in the table.
(b) **The row labels are the VARIABLE HELD or SWEPT, not the sample.** Under "Autoclaved", the rows "29 % H2O"
and "40 % H2O" are (K at 121 C − K at 100 C) at that moisture; the rows "100 °C" and "121 °C" are (K at 40 % −
K at 29 %) at that temperature. So the **"40 % H2O" row being positive means raising temperature from 100 to
121 C RAISES K when the protein is wet**, and the **"100 °C" row being negative means raising moisture from
29 % to 40 % LOWERS K when the treatment is at 100 C**. The paper states both in words on p. 105: an increase
in temperature at 40 % moisture increases K by **56-59 % for ketones, 64-81 % for alcohols, 64 % for hexanal**;
an increase in moisture at 100 C decreases K by **33-47 %**.

### Table 7 (p. 111). "Comparison of heats of adsorption on three solid substrates. The present values are averages for all treatments of Edi-Pro A."

Values are **−ΔH in kcal/mol^1**.

| Flavor | Aspelund, 1978 — Edi-Pro A | Current Study — Heat-Treated Edi-Pro A | McMullin et al., 1975 — Lactose |
|---|---:|---:|---:|
| Pentanone | 2.92 `[C]` | 8.10 `[F]` | 10.81 `[C]` |
| Hexanone | 6.04 `[C]` | 9.81 `[F]` | 12.03 `[C]` |
| Heptanone | 8.11 `[C]` | 11.56 `[F]` | 13.09 `[C]` |
| Hexanal | 8.89 `[C]` | 9.90 `[F]` | 11.66 `[C]` |
| Pentanol | 11.25 `[C]` | 14.50 `[F]` | 15.70 `[C]` |
| Hexanol | 13.89 `[C]` | 16.45 `[F]` | 17.57 `[C]` |
| Heptanol | 18.06 `[C]` | **17.80** `[F]` | 17.18 `[C]` |

Footnote exactly as printed: *^1 All values for each work were averaged.*

**The Aspelund column is the same data as `aspelund1983_extraction.md` Table I, from the 1978 thesis.** Six of
the seven values — 6.04, 8.11, 8.89, 11.25, 13.89, 18.06 — are **identical to Aspelund 1983's Table I to the
last digit**. The seventh, pentanone at **2.92**, is a *measured* thesis value for a compound that Aspelund
1983 declared statistically non-significant and reported only as an **extrapolated 4.3** (Table II, footnote a).
So this table preserves a measured number the 1983 paper suppressed (Flags 8).

**Internal check (mine): the "Current Study" column against the six rows of Table 4.** Averaging the four
autoclaved rows plus untreated plus sheared: pentanone **8.102** vs 8.10 ✓; hexanone **9.807** vs 9.81 ✓;
heptanone **11.56** vs 11.56 ✓; hexanal **9.887** vs 9.90 ✓; pentanol **14.50** vs 14.50 ✓; hexanol **16.45**
vs 16.45 ✓; **heptanol 17.00 vs a printed 17.80 — a 0.80 kcal/mol MISMATCH, the only one in the table.** Six
of seven reproduce exactly, so the averaging basis is certain and heptanol is a genuine error somewhere. To
reach 17.80 the six heptanol values would have to sum to 106.8 rather than the 102.0 they do. **Heptanol is
also the one compound whose ΔH rests on only two temperatures in three of the runs** (Table 4, footnote 1).
**Treat the heptanol ΔH column as unreliable in both tables.**

### Numbers printed in the running text

| quantity | value | where | class |
|---|---|---|---|
| Edi-Pro A composition (manufacturer) | **93.5 % protein, 5.5 % moisture, 0.2 % fat** | p. 100 | `[C]` (Ralston-Purina, not assayed here) |
| autoclave hold | **20 min** in a preheated autoclave, protein spread thinly on trays | p. 101 | `[M]` |
| moisture conditioning | water added, dried back at **55 C** in an air-circulation drier with occasional stirring | p. 101 | `[M]` |
| extruder | Wenger X-5, **700 rpm, no die**, feed at **33 % moisture**, **maximum barrel temperature 71 C** | p. 101 | `[M]` |
| packing sizing | **230/325 mesh**, Tyler sieves | p. 101 | `[M]` |
| packing mass | **~1.4 g treated, 0.5 g unprocessed**, at 50 psig head pressure | p. 101 | `[M]` |
| n-hexane | **"not significantly retained"** — excluded from all tables | p. 103 | `[M]` (a measured null) |
| K at 60 and 80 C | **"not shown"**; "differ in magnitude but qualitatively show the same effects of processing"; complete data in Crowther 1979 (thesis) | p. 103 | **not available** |
| peak tailing | minimised by reducing injected amount but **"cannot be eliminated for alcohols"** | p. 103 | `[M]` |
| **effect of raising 100 → 121 C at 40 % moisture** | K **increases**: **56-59 % (ketones), 64-81 % (alcohols), 64 % (hexanal)** | p. 105 | `[M]` |
| **effect of raising 29 → 40 % moisture at 100 C** | K **decreases 33-47 %** | p. 105 | `[M]` |
| **effect of autoclaving vs untreated** | K falls **47-49 % (ketones), 35-40 % (alcohols), 43 % (hexanal)** | p. 111 | **`[M]` — the headline** |
| **effect of shear vs untreated** | significant **only on the alcohols**, values **increased 14-46 %** | p. 111 | `[M]` |
| lowest K of all | at **121 C, 29 % H2O** | p. 111 | `[M]` |
| alcohol − carbonyl ΔH gap | **6.2-6.6 kcal/mol** | p. 111 | `[M]` |
| n-heptane ΔH as a van der Waals reference | **6.38 kcal/mol** on lactose | p. 111 | `[C]` (McMullin 1975) |
| the same, from Aspelund | **6.48** | p. 112 | `[C]` (Aspelund 1978) |
| heptane→heptanone and heptanone→heptanol ΔH steps | **5-6 kcal/mol** each | p. 112 | `[C]`/interpretation |
| binding model proposed | non-specific interaction **plus one hydrogen bond for carbonyls, two for alcohols**; requires the protein to have "numerous polar binding sites" | p. 112 | interpretation |
| **the mechanism claim** | *"Since ΔH values were unaffected by processing, the observed changes in K are postulated to be the result of changes in the NUMBER OF BINDING SITES"* | p. 112 | **`[F]`/interpretation — the sentence the repository cares about** |
| **the dry-vs-aqueous reversal** | *"The comparative extents of binding (as indicated by K) in this dry system are predominantly the reverse of the order found in aqueous solutions of soy proteins"* | p. 112 | **interpretation — a hard transfer prohibition** |
| aqueous ordering, for contrast | binding decreases **aldehyde > ketone > alcohol** | p. 100 | `[C]` (Beyeler & Solms 1974) |
| aqueous soy, alcohols do not bind at all | Gremli: alcohols did not bind; aldehydes bind **both reversibly and irreversibly**; ketones bind reversibly | p. 100 | `[C]` (Gremli 1974, 5 % soy protein solution) |
| **reversible-only property of the method** | *"the only molecules detected were those that traversed the column and thus were only physically (reversibly) adsorbed"* | p. 112 | **`[M]` by construction** |
| Arai's contrary result | native / partly denatured / denatured soy retained **increasing** amounts of hexanal and hexanol, and **more hexanal than hexanol** in every case — attributed by Crowther to chemical (irreversible) bonding | p. 112 | `[C]` (Arai 1970) |
| the denaturation picture | **chain unfolding** exposing nonpolar regions, decreasing solubility and the availability of polar hydrogen-bonding sites; heat aggregation with intermolecular peptide bonds also possible | pp. 112-113 | interpretation (Solms 1973, Burgess & Stanley 1976) |
| SEM verdict | **"no evident visual difference between the sheared and heat-treated protein"** (Fig. 4: unprocessed 930×, heat-treated before grinding 225×, sheared 290×, heat-treated 700×) | p. 112 | `[M]` (morphological only) |
| **the variable the authors say they should have measured** | *"This suggests that we ought to consider WATER ACTIVITY as the important process variable. However, such measurements were not made."* | p. 114 | — |
| correlation coefficients of the ΔH regressions (Fig. 3, carbonyls, one 121 C/29 % packing) | heptanone **0.9930**, hexanone **0.9973**, hexanal **0.9978**, pentanone **0.9378** | Fig. 3, p. 110 | `[M]` (printed inside the panel) |

**Figure-only quantities.** The Fig. 2 alcohol correlation coefficients are inside the panel but are not
legible in this scan and are **not typed here**. Fig. 1's nine peaks carry no numbers ("The order of the peaks
has no significance"). Fig. 4's micrographs carry only magnifications. The 60 C and 80 C adsorption
coefficients exist and are **not printed anywhere in this paper**.

### Arithmetic on the printed constants (all mine)

**1. The heat-treatment effect as a ratio, which is the form the registry uses.** Untreated K divided by the
mean of the four autoclaved K's, at 70 C:

| compound | untreated K | autoclaved mean K | **untreated / autoclaved (mine)** | % fall (mine) | paper's stated % (p. 111) |
|---|---:|---:|---:|---:|---|
| 2-pentanone | 1.56 | 0.825 | **1.89x** | 47.1 % | 47-49 % (ketones) |
| 2-hexanone | 3.48 | 1.78 | **1.96x** | 48.9 % | " |
| 2-heptanone | 6.73 | 3.43 | **1.96x** | 49.0 % | " |
| 1-pentanol | 11.3 | 7.40 | **1.53x** | 34.5 % | 35-40 % (alcohols) |
| 1-hexanol | 21.7 | 14.1 | **1.54x** | 35.0 % | " |
| 1-heptanol | 52.1 | 31.3 | **1.67x** | 39.9 % | " |
| hexanal | 2.77 | 1.58 | **1.75x** | 43.0 % | 43 % |

**Every recomputed value falls inside the paper's own stated bands.** The ketone class is the most affected and
the alcohol class the least, with the aldehyde between — a **class-dependent** preheat factor, not a single
number.

**2. The shear effect as a ratio (mine).** Sheared K over untreated K at 70 C: pentanone **1.045x**, hexanone
**0.931x**, heptanone **0.966x**, pentanol **1.460x**, hexanol **1.369x**, heptanol **1.140x**, hexanal
**1.090x**. The three alcohols reproduce the paper's "increased from 14-46 %" exactly. **Two of the three
ketones went DOWN**, and hexanone's fall (−6.9 %) is one of the starred entries in Table 6.

**3. Within the autoclave block, the four conditions span nearly 2x (mine).** Highest over lowest K among the
four autoclave conditions: pentanone **1.88x**, hexanone **1.78x**, heptanone **1.50x**, pentanol **1.64x**,
hexanol **1.77x**, heptanol **1.81x**, hexanal **1.74x**. **So "autoclaved" is not one state**: 121 C at 40 %
moisture retains nearly twice the flavour of 100 C at 40 % moisture. The lowest K throughout is **121 C, 29 %**
except for pentanone, hexanone, heptanone, pentanol, hexanol, heptanol and hexanal where Table 3 makes **100 C,
40 %** the lowest — **the paper's own claim on p. 111 that "the lowest K values occurred at 121 °C, 29 % H2O"
is CONTRADICTED by its own Table 3 for all seven compounds** (Flags 1).

**4. The chain-length slope, per treatment (mine).** Consecutive K ratios within a homologous series at 70 C:

| series and treatment | consecutive ratios | geometric mean per CH2 |
|---|---|---:|
| 2-alkanones, untreated (C5→C6→C7) | 2.231, 1.934 | **2.077** |
| 2-alkanones, autoclaved mean | 2.158, 1.927 | **2.039** |
| 2-alkanones, sheared | 1.988, 2.006 | **1.997** |
| 1-alkanols, untreated (C5→C6→C7) | 1.920, 2.401 | **2.147** |
| 1-alkanols, autoclaved mean | 1.905, 2.220 | **2.057** |
| 1-alkanols, sheared | 1.800, 2.000 | **1.897** |

**Two findings.** (a) The slope is **1.90-2.15x per CH2** on soy at 70 C, against the shipped
`CHAIN_LENGTH_SLOPE_PER_CH2 = 2.81` — **24-32 % low**, in the same direction as Aspelund's 2.23-2.27x at 90 C.
(b) **The slope is essentially INVARIANT to the treatment**: 2.08 → 2.04 → 2.00 for the ketones across
untreated / autoclaved / sheared, and 2.15 → 2.06 → 1.90 for the alcohols. **A preheat moves the level and
leaves the chain-length structure alone.** That is exactly the shape a "number of sites falls, strength per site
does not" model predicts, and it is an independent confirmation of the paper's ΔH null through a different
quantity.

**5. Per-gram retention, because K is per SQUARE METRE and the registry is per GRAM (mine).** K = V_m / S, so
V_m = K × S in mL/g. Using Table 2's surface areas (autoclaved conditions averaged over their two replicate
determinations: 0.175, 0.22, 0.21, 0.24 m^2/g; untreated 0.20; sheared 0.28):

| compound | untreated V_m | autoclaved-mean V_m | **untreated / autoclaved, per GRAM (mine)** | sheared V_m | **sheared / untreated, per GRAM (mine)** |
|---|---:|---:|---:|---:|---:|
| 2-pentanone | 0.312 mL/g | 0.173 | **1.81x** | 0.456 | **1.46x** |
| 2-hexanone | 0.696 | 0.372 | **1.87x** | 0.907 | **1.30x** |
| 2-heptanone | 1.346 | 0.724 | **1.86x** | 1.820 | **1.35x** |
| 1-pentanol | 2.260 | 1.560 | **1.45x** | 4.620 | **2.04x** |
| 1-hexanol | 4.340 | 2.957 | **1.47x** | 8.316 | **1.92x** |
| 1-heptanol | 10.420 | 6.567 | **1.59x** | 16.632 | **1.60x** |
| hexanal | 0.554 | 0.331 | **1.67x** | 0.846 | **1.53x** |

**The autoclave conclusion survives the change of basis almost unchanged** (1.45-1.87x per gram against
1.53-1.96x per m^2) — because the surface areas happen to be similar. **The shear conclusion does NOT.** Shear
raised the surface area from 0.20 to 0.28 m^2/g (+40 %), so on a per-gram basis shear increases retention for
**every compound including the ketones** — hexanone goes from **0.93x per m^2 to 1.30x per gram**, a reversal
of sign. **The registry works per gram. If any Crowther ratio is ever carried, it must be the per-gram one, and
the shear rows change direction between the two bases** (Flags 9).

**6. This paper against its own predecessor, and they do not agree (mine).** Aspelund 1983 measured V_S on
**untreated Edi-Pro A** in **the same units (mL/m^2)** by **the same construction** at **90 C**. Crowther
measured K on untreated Edi-Pro A at **70 C**. Adsorption must fall with temperature, so Aspelund's 90 C values
should be **below** Crowther's 70 C values. Projecting Crowther's untreated K from 70 C to 90 C using
**Crowther's own untreated ΔH** and his own Eq. 4 form (ln(V_m/T) = −ΔH/RT + A):

| compound | Crowther K at 70 C | Crowther ΔH untreated | **projected K at 90 C (mine)** | Aspelund V_S at 90 C | **Aspelund / projected (mine)** |
|---|---:|---:|---:|---:|---:|
| 2-hexanone | 3.48 | −11.0 | 1.52 | 3.29 | **2.2x** |
| 2-heptanone | 6.73 | −12.7 | 2.55 | 7.62 | **3.0x** |
| hexanal | 2.77 | −11.1 | 1.20 | 3.14 | **2.6x** |
| 1-pentanol | 11.3 | −13.1 | 4.15 | 10.99 | **2.6x** |
| 1-hexanol | 21.7 | −16.7 | 5.96 | 31.24 | **5.2x** |
| 1-heptanol | 52.1 | −16.8 | 14.19 | 70.95 | **5.0x** |

**Same laboratory, same commercial product, same method family, same units — and the two determinations
disagree by 2.2 to 5.2x**, in the direction of Aspelund reading systematically high. Crowther's corrections go
the wrong way to close the gap: he applies a **James-Martin compressibility factor j** (less than 1) and a
column/flowmeter temperature ratio that Aspelund does not, which would make Crowther's K *smaller* still
relative to an uncorrected Aspelund V_S. And the two papers also disagree on the **slope**: Crowther's untreated
ΔH values are **1.4-1.8x more negative** than Aspelund's for the ketones and hexanal (11.0 vs 6.04, 12.7 vs
8.11, 11.1 vs 8.89), which Crowther notes ("Values calculated in our work are generally higher (i.e. more
negative) than Aspelund's", p. 111) without explaining. **Two papers from one laboratory on one product differ
by up to 5x on the level and up to 1.8x on the temperature slope. Neither absolute scale should be trusted;
only within-study ratios should be carried.** (This is the same conclusion Amendment 4 reached about Meynier's
and Leksrisompong's static-headspace scales, reached here by a different route.)

**7. The functional-group contrasts at six carbons, untreated, 70 C (mine).** hexanol / hexanal = **7.83x**;
hexanol / hexanone = **6.24x**; hexanone / hexanal = **1.26x**. Aspelund's dry-phase equivalents at 90 C were
9.9x, 9.5x and 0.95x. **The two studies agree that the alcohol dominates by roughly an order of magnitude and
that the ketone and the aldehyde are close; they disagree on which of the two carbonyls is larger.** Both
orderings are the **reverse** of the aqueous ordering the same paper cites (p. 112). None of these is
transferable.

**8. What Table 6's differences look like as the ratios the paper reports in words (mine).** For "40 % H2O"
(the effect of 100 → 121 C at 40 % moisture): pentanone 0.879/0.556 = **1.58x**, hexanone 1.91/1.20 = **1.59x**,
heptanone 3.94/2.62 = **1.50x** (the paper says 56-59 % for ketones — note heptanone at 50 % sits just outside
the stated band), pentanol 9.04/5.51 = **1.64x**, hexanol 17.0/9.6 = **1.77x**, heptanol 38.4/21.2 = **1.81x**
(paper: 64-81 %), hexanal 1.79/1.09 = **1.64x** (paper: 64 %). **Six of seven reproduce the paper's stated
percentages exactly; heptanone is 50 % against a stated 56-59 %**, which is a rounding-band slip in the text
rather than a table error.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** `hexanal` is keyed. `1_hexanol` is keyed. `2_heptanone` is
keyed. **`2_pentanone`, `2_hexanone`, `1_pentanol` and `1_heptanol` are not**, and none of them is a panel
target; nothing here requires a new key, because no row is shippable as a constant. `COMPOUND_STRUCTURE` in
`parameters_matrix.py` carries `hexanal` and `2_heptanone` and would need an `alcohol` class extension for the
1-alkanols, which Wave B26 began with `z_2_penten_1_ol`.

Every row below shares: **Edi-Pro A soy protein isolate (93.5 % protein per the manufacturer), processed as
stated, then dried back to the original moisture, ground, sieved 230/325 mesh, packed DRY as a GC stationary
phase; nitrogen carrier at 23-28 mL/min, 50 psig head; K measured at 70 C; ΔH from 60/70/80 C; inverse gas
chromatography, reversible adsorption only.** There is **no solvent, no pH and no protein loading in g/L** on
any of them.

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **adsorption coefficient K, hexanal, UNTREATED soy** | **2.77 ± 0.45** | mL/m^2 | 70 C, dry, inverse GC | Table 3, p. 106 | `binding_constant` (gas-solid; NOT an aqueous K) |
| adsorption coefficient K, hexanal, autoclaved 121/29, 121/40, 100/29, 100/40 | **1.54 ± 0.11 / 1.79 ± 0.27 / 1.90 ± 0.30 / 1.09 ± 0.02** | mL/m^2 | as above, after 20 min autoclave | Table 3, p. 106 | `binding_constant` (gas-solid) |
| adsorption coefficient K, hexanal, sheared | **3.02 ± 0.07** | mL/m^2 | as above, extruded 71 C | Table 3, p. 106 | `binding_constant` (gas-solid) |
| adsorption coefficient K, 2-heptanone, untreated / autoclaved-mean / sheared | **6.73 ± 0.53 / 3.43 ± 0.32 / 6.50 ± 0.26** | mL/m^2 | 70 C | Table 3, p. 106 | `binding_constant` (gas-solid) |
| adsorption coefficient K, all 7 compounds × 6 treatments | see Table 3 | mL/m^2 | 70 C | Table 3, p. 106 | `binding_constant` (gas-solid) |
| **PREHEAT FACTOR on K, ketones (autoclave vs untreated)** | **0.51-0.53** (a 47-49 % fall) | × | 100-121 C, 20 min, 29-40 % moisture | p. 111 and Table 3 (mine) | **`within_study_ratio`** — the number the `matrix_sites` gap wants |
| **PREHEAT FACTOR on K, alcohols** | **0.60-0.66** (a 35-40 % fall) | × | as above | p. 111 and Table 3 (mine) | **`within_study_ratio`** |
| **PREHEAT FACTOR on K, hexanal** | **0.57** (a 43 % fall) | × | as above | p. 111 and Table 3 (mine) | **`within_study_ratio`** |
| the same three, **on a PER-GRAM basis** | ketones **0.53-0.55**, alcohols **0.63-0.69**, hexanal **0.60** | × | as above | K × S from Tables 2 and 3 (mine) | `derived_assumption` — the basis the registry actually uses (Flags 9) |
| **SHEAR FACTOR on K, alcohols** | **1.14-1.46** (a 14-46 % RISE) | × | extruder, 700 rpm, 33 % moisture, max 71 C | p. 111 and Table 3 (mine) | **`within_study_ratio`** — a processing effect with the OPPOSITE sign |
| shear factor on K, ketones and hexanal | **0.93-1.09** | × | as above | Table 3 (mine) | `within_study_ratio` (two of four below 1) |
| shear factor **per gram** | **1.30-2.04, every compound above 1** | × | as above | K × S (mine) | `derived_assumption` — **sign flips vs the per-m^2 basis** (Flags 9) |
| effect of 100 → 121 C at 40 % moisture | K rises **1.50-1.81x** (paper: 56-59 % ketones, 64-81 % alcohols, 64 % hexanal) | × | 20 min autoclave | p. 105, Table 6 p. 109 (mine) | `within_study_ratio` |
| effect of 29 → 40 % moisture at 100 C | K falls **33-47 %** | % | 20 min autoclave | p. 105, Table 6 p. 109 | `within_study_ratio` |
| spread across the four autoclave conditions | **1.50-1.88x** high over low | × | — | Table 3 (mine) | `within_study_ratio` — "autoclaved" is not one state |
| **heat of adsorption ΔH is UNCHANGED by processing** | ANOVA p = **0.076 to 0.974**, nothing below 0.05 | — | 2 × 2 factorial, n = 4 conditions | Table 5, p. 108 | **`measured_bound`** — a NULL result, and the paper's central claim (Flags 2) |
| heat of adsorption ΔH, untreated soy, 7 compounds | **−10.1 / −11.0 / −12.7 / −13.1 / −16.7 / −16.8 / −11.1** (pentanone, hexanone, heptanone, pentanol, hexanol, heptanol, hexanal) | kcal/mol | van 't Hoff over 60-80 C | Table 4, p. 107 | `binding_constant` — **THERMODYNAMIC, NOT AN ACTIVATION ENERGY** (Flags 10) |
| heat of adsorption ΔH, all 7 compounds × 6 treatments | see Table 4 | kcal/mol | as above | Table 4, p. 107 | `binding_constant` — thermodynamic |
| chain-length slope, 2-alkanones, untreated / autoclaved / sheared | **2.08 / 2.04 / 2.00** | × per CH2 | 70 C, dry | Table 3 (mine) | **`within_study_ratio`** — vs shipped 2.81 |
| chain-length slope, 1-alkanols, untreated / autoclaved / sheared | **2.15 / 2.06 / 1.90** | × per CH2 | as above | Table 3 (mine) | **`within_study_ratio`** |
| **the chain-length slope is invariant to preheat** | moves by ≤ 4 % (ketones) and ≤ 12 % (alcohols) across treatments | — | as above | Table 3 (mine) | **`within_study_ratio`** — a useful null |
| hexanol / hexanal and hexanol / hexanone contrasts, untreated | **7.83x / 6.24x** | × | 70 C, dry | Table 3 (mine) | `within_study_ratio` — **dry-phase ordering, the REVERSE of aqueous** (Flags 8 of Aspelund; here p. 112) |
| **protein dispersibility index, neutralised** | **88.4 (untreated) / 48.7 (sheared) / 7.8 (autoclaved, all four pooled)** | % | — | Table 2, p. 105 | **`level_only`** — the denaturation gradient that pairs with the K fall |
| protein dispersibility index, isoelectric | 1.3 / 1.6 / 1.9 | % | — | Table 2, p. 105 | `level_only` |
| disc gel electrophoresis | several strong bands / less distinct with some missing / **no bands** | — | untreated / sheared / autoclaved | Table 2, p. 105 | `level_only` |
| specific surface area | untreated **0.20**, sheared **0.28**, autoclaved **0.17-0.28** | m^2/g | argon single-point BET ±3 % (one krypton multipoint ±1 %) | Table 2, p. 105 | `measured_bound` |
| Edi-Pro A protein content | **93.5 %** (5.5 % moisture, 0.2 % fat) | % | manufacturer's figure | p. 100 | `level_only` — **`[C]`, not assayed here** |
| n-hexane retention on soy | **not significantly retained** | — | 70 C, dry | p. 103 | **`measured_bound`** — a measured null; alkanes do not adsorb detectably |
| per-gram retention volumes V_m = K × S | see section 3 item 5 | mL/g | 70 C, dry | Tables 2 and 3 (mine) | **`derived_assumption` — DO NOT SHIP.** Same units family as K_g, different quantity |
| cross-paper disagreement with Aspelund on the same product | **2.2 to 5.2x** on the level, **1.4-1.8x** on the ΔH slope | × | same lab, same Edi-Pro A, same units | section 3 item 6 (mine) | **`within_study_ratio`** — the strongest same-lab scale warning in the dry-phase literature |

### Can these be put on the same basis as the shipped binding constants, i.e. converted to K_g in L/g?

**No, for the same structural reason as Aspelund 1983, and one more besides.**

- **There is no water leg.** `K_g = (K_water/K_matrix − 1) / protein_g_per_L` needs an air/water coefficient
  and an air/matrix coefficient from the same system. This experiment has nitrogen as the mobile phase and dry
  protein as the stationary phase. There is nothing to subtract, so the *excess-over-water* that K_g measures
  cannot be formed. That is what inverse GC is, not an omission.
- **There is no protein loading in g/L.** The paper prints a packing mass (1.4 g treated, 0.5 g untreated) and
  a surface area (0.17-0.28 m^2/g). Neither is a concentration. Section 3 item 5 converts K to a per-gram
  retention volume anyway, and the result is a number in mL/g that is **not** K_g.
- **And the absolute scale is independently untrustworthy**, by this paper's own comparison with its
  predecessor: 2.2-5.2x disagreement between two runs of the same method on the same product in the same
  laboratory (section 3 item 6). Even if a conversion existed, the absolute value should not be shipped.
- **A molar mass is not needed.** K is already per gram of packing once multiplied by S.

**What IS on a common basis, and it is the point of reading this paper: the RATIOS.** The registry's whole
`REVERSIBLE_BINDING` construction under Amendment 4 is "carry the within-study ratio, refuse the absolute
scale", because the absolute static-headspace scale was 6.24x low for Meynier and 6-17x low for Leksrisompong.
Crowther's absolute scale is suspect by 2.2-5.2x against its own predecessor and his ratios are within-run:
same column, same instrument, same day, same compound, treated against untreated. **The preheat factor
(0.51-0.66 on K, 0.53-0.69 per gram) and the shear factor (1.14-1.46 on alcohols) are the transferable
objects.** They are not binding constants and they should not be shipped as one; they are **a candidate
multiplier on a site count**, which is what `matrix_sites.py` charges.

**What a defensible use would look like, and what it would need.** `matrix_sites.py` charges free thiol,
disulfide and amine densities in mmol/g at the start of the cook. Crowther says a 20-min autoclave at 100-121 C
removes 35-49 % of the *flavour-adsorption* sites on soy while leaving their individual strength unchanged. To
use that as a preheat multiplier on the covalent site pools would require **three things this paper does not
supply**: (i) that the polar hydrogen-bonding sites an inverse-GC K counts are the same sites a lysine-amine
adduct channel consumes — they are almost certainly not, since Crowther's mechanism is a *hydrogen bond* and
`matrix_sites`' channel is a *covalent Schiff-base/Michael* addition; (ii) an aqueous replication; (iii) a
dependence on time-at-temperature, since only a single 20-min hold was tested. **The honest use is as a bound
and a direction: preheating soy at process temperature reduces its reversible flavour retention by roughly
half, class-dependent, and does so by removing sites rather than weakening them.** That is a `measured_bound`
plus a `within_study_ratio`, not a parameter.

**Nothing here goes to `matrix_sites.py` as a rate.** No rate constant, no time series, no activation energy
appears in this paper. See Flags 10.

## 5. Flags

1. **The paper contradicts its own Table 3 about where the minimum is.** p. 111 states "The lowest K values
   occurred at 121 °C, 29 % H2O (see Table 3 for K values)". **Table 3 says otherwise for all seven compounds**:
   the 100 C / 40 % row is lower than the 121 C / 29 % row in every column (pentanone 0.556 vs 0.822; hexanone
   1.20 vs 1.89; heptanone 2.62 vs 3.49; pentanol 5.51 vs 6.66; hexanol 9.6 vs 13.5; heptanol 21.2 vs 30.4;
   hexanal 1.09 vs 1.54). The table is internally consistent and its Mean row closes, so **the sentence is
   wrong, not the table**. Anyone reading the discussion without the table will take away the wrong worst case.
2. **The central claim — "ΔH unchanged" — is a NULL RESULT on four conditions with large error bars.** Table 5
   reports no ΔH p-value below 0.076, and the t-tests were never run on ΔH at all "since the data clearly
   indicated that there was no significant change" (p. 103). But Table 4's standard errors reach **±3.20**
   (hexanone at 100 C / 40 %) and **±1.78** (heptanol at 121 C / 29 %), and the **untreated and sheared rows
   carry no error bars whatsoever**, so the comparison that matters most — treated against untreated — has no
   variance estimate on one side. Meanwhile the *point estimates* do move: the pentanone ΔH ranges from −6.26
   to −10.0 across the four autoclave conditions, a 60 % spread. **"Not significant" here means "not resolvable
   at n = 4", not "measured to be unchanged".** The inference that heating changes site count and not site
   strength rests on this null and should be reported with its power, not as a fact.
3. **Table 6's significance stars run BACKWARDS from convention.** *"\*P < 0.01; \*\* P < 0.05; \*\*\* P < 0.1;
   \*\*\*\* P < 0.2"* — one star is strongest, four weakest. A reader who assumes the usual convention will
   invert the entire table's evidential weight. Note also that the quantity being starred is a **difference**
   (ΔK), not a ratio, so the stars are not comparable across compounds of different magnitude: hexanal's
   −1.19\*\*\* is a 43 % effect and heptanol's −20.8\*\* is a 40 % effect, but the raw numbers differ 17-fold.
4. **This paper and Aspelund 1983 disagree by 2.2-5.2x on the same product with the same method.** Detailed in
   section 3 item 6. Both are Iowa State, both use Edi-Pro A, both report a per-square-metre gas-solid constant,
   and once the 20 C temperature difference is removed using Crowther's own heats of adsorption the two
   determinations differ by a factor of two to five, with Aspelund high. The **ordinate of the van 't Hoff
   regression also differs** — Aspelund regresses ln(t_cor) against 1/T (p. 540 of that paper), Crowther
   regresses ln(V_m/T) against 1/T (Eq. 4, p. 102) — which changes the fitted ΔH by roughly R·T ≈ 0.7 kcal/mol
   and is nowhere near enough to explain the 1.4-1.8x gap in ΔH. **Neither paper's absolute scale is usable.
   Both papers' internal ratios remain usable.**
5. **The untreated control was packed at a third of the mass of the treated samples.** "Approximately 1.4 g of
   the treated samples and **0.5 g of unprocessed protein** could be used while attaining reasonable flow rates"
   (p. 101). K is normalised by mass and by surface area, so in principle this cancels; in practice a column
   with 0.5 g of packing has a different void fraction, a different pressure drop (and therefore a different
   James-Martin j), and a different residence-time distribution from one with 1.4 g. **The single most important
   comparison in the paper — treated against untreated — is also the one comparison in which the columns were
   not physically alike.** The paper does not test for a mass artefact; Aspelund 1983, on a similar question,
   ran an explicit mass-temperature interaction test and found none, but that was over a 1.4-1.5 g range, not a
   0.5-1.4 g range.
6. **The injected quantity is never stated and was adjusted per compound by an undocumented procedure.**
   "Where more dilute samples were required, the syringes were flushed several times before injection" (p. 101).
   Aspelund 1983 quantified its injections at 10^-9 to 10^-6 g; this paper does not. The linear-isotherm
   (Henry's-law) assumption on which Eq. 4 and the definition of K both depend is asserted — "Assuming that we
   are working in the linear portion of the adsorption isotherm (very low gas concentrations)" — and never
   tested. And the paper concedes that **peak tailing "cannot be eliminated for alcohols"** (p. 103), which is
   the classic signature of a *non-linear* isotherm on a heterogeneous surface. **The alcohols are the compound
   class with (a) tailing peaks, (b) the largest K's, (c) the anomalous shear response and (d) the heptanol
   arithmetic error in Table 7. Treat the alcohol column as the weakest in the paper.**
7. **The method cannot see irreversible binding, which is a strength and a blind spot.** "The only molecules
   detected were those that traversed the column and thus were only physically (reversibly) adsorbed" (p. 112).
   So K is a clean reversible constant — unlike a headspace-depletion measurement, which counts covalent capture
   as binding (the reason Bi 2022's (E)-2-octenal row is quarantined in `REVERSIBLE_BINDING`). But it also means
   **the paper is blind to precisely the channel `matrix_sites.py` models.** If autoclaving *created* covalent
   aldehyde-binding capacity while destroying hydrogen-bonding capacity, this experiment would report only the
   fall. The paper itself raises the possibility against Arai 1970, which found *more* retention on *more*
   denatured soy and which Crowther attributes to chemical bonding (p. 112). **Note the registry's
   `kg_hexanal_soy_denatured` (1.47e-3 L/g) traces to Arai 1970 via Damodaran — i.e. to the study Crowther
   argues was measuring a different phenomenon from his own. The two are not alternative estimates of one
   quantity and must not be averaged.**
8. **"Autoclaved" pools four materially different states, and the PDI is pooled too.** The four conditions span
   1.5-1.9x in K (section 3 item 3), and their surface areas span 0.17-0.28 m^2/g. Yet Table 2 reports a
   **single** PDI (1.9 isoelectric, 7.8 neutralised) and a single electrophoresis result ("No bands") for all
   four. So the paper cannot correlate degree of denaturation with degree of binding loss **within** the
   autoclave block — the very correlation its discussion asserts ("both electrophoresis and PDI results changed
   most drastically for the moist-heat treated samples ... This, in turn, coincides with the reduction in
   adsorption coefficients", p. 113). **The coincidence is asserted across three coarse groups, not measured
   across a gradient.**
9. **Per-square-metre and per-gram give different answers for shear, including a sign change.** K is normalised
   by surface area by design ("independent of packing mass or surface area", p. 103), and shear raised the
   surface area 40 % (0.20 → 0.28 m^2/g). On K, sheared hexanone is **0.93x** untreated; per gram it is
   **1.30x**. **The registry stores per-gram constants. Whichever basis is used must be stated, and the shear
   result is basis-dependent.** The autoclave result is not — it survives the change of basis (1.45-1.87x
   per gram against 1.53-1.96x per m^2) because the autoclaved surface areas happen to bracket the untreated one.
10. **The heats of adsorption are van 't Hoff quantities, not activation energies.** Eq. 4 (p. 102) fits ΔH from
    the slope of **ln(V_m/T) against 1/T**, an equilibrium retention volume against reciprocal temperature,
    after Kiselev & Yashin 1969. These are **isosteric heats of a reversible physical adsorption**, printed with
    a negative sign because adsorption is exothermic. **They are NOT E_a**, they must not enter
    `matrix_sites.py`'s `ea_band_kj_mol`, and their magnitudes invite the mistake: −11.1 kcal/mol for hexanal is
    **−46.4 kJ/mol (mine)**, which would look like a plausible activation energy and is not one. An activation
    energy cannot be negative. **There is no rate and no time axis anywhere in this paper**, only a 20-minute
    treatment hold that is never varied.
11. **A single time-at-temperature was tested.** Every autoclave treatment is **20 min**. There is no 5-minute
    or 60-minute arm, so nothing here says whether the 35-49 % loss is an early plateau or the midpoint of a
    continuing decline. A model that needs to know what a 3-minute or a 40-minute hold does gets no guidance.
12. **Water activity is the variable the authors themselves say should have been measured and was not.**
    p. 114: *"This suggests that we ought to consider water activity as the important process variable. However,
    such measurements were not made."* The moisture-temperature **interaction** is the only factor that reached
    significance on K (Table 5: 0.025 and 0.047 for pentanone and hexanone), and it is precisely the interaction
    an a_w treatment would collapse into one variable. **The design's dominant effect is one the design cannot
    interpret.**
13. **The shear arm is not a temperature-matched control.** The extruder reached a maximum barrel temperature of
    **71 C** at 33 % moisture, against autoclave holds at 100-121 C at 29-40 %. So "shear" differs from
    "autoclave" in temperature, moisture, mechanical energy and residence time simultaneously. **The +14-46 %
    alcohol response cannot be attributed to shear as opposed to mild heat**, and its PDI (48.7 neutralised)
    sits squarely between untreated (88.4) and autoclaved (7.8) — consistent with it being simply a milder
    denaturation, not a different kind.
14. **Table 7's heptanol value does not reproduce from Table 4** — printed 17.80, recomputed 17.00 (section 3).
    Six of the seven rows reproduce exactly, so the averaging basis is certain. Heptanol is also the compound
    with only two temperatures in three of its runs (Table 4, footnote 1). **The heptanol ΔH column should be
    treated as unreliable in both tables**, and the discrepancy noted if the paper is cited for it.
15. **Two internal typographical inconsistencies to note when citing.** The Table 1 caption says "Varian 3740"
    while the Methods say "Varian 3700" (pp. 101-102). Table 4's column header prints "Hepanol" for heptanol.
    Neither affects a number. Separately, the **publication year is ambiguous** — journal footer 1980, copyright
    1981, and the same laboratory's own 1983 paper cites it as 1981 (section 0).
16. **The isolate is a 1980s commercial product characterised only by the manufacturer's own composition
    figure.** No thiol, disulfide or free-amine assay; no molar mass; no independent protein determination.
    `data/species/protein_matrices.yml` charges soy site densities from other, modern preparations. Pairing a
    Crowther preheat ratio with those densities would be a cross-preparation pairing 45 years apart and should be
    labelled as one.
17. **No DOI is printed** (section 0). Cite by volume/page and note the year ambiguity.
18. **What this paper does NOT contain**: any aqueous phase, pH, ionic strength or buffer; any binding constant
    in M^-1 or L/g; any protein concentration in g/L; any rate constant or activation energy; any covalent adduct
    measurement; any 2-alkenal, pyrazine, pyridine, furan or sulfur compound; any water-activity measurement; any
    time-at-temperature series; **the 60 C and 80 C adsorption coefficients (stated to exist, referred to an
    unpublished 1979 thesis, not printed here)**; any standard error on the untreated or sheared heats of
    adsorption; any sensory measurement; any odour threshold; any supplementary material.
19. **What to request** (Iowa State, 1980 — recorded for completeness): (i) the 60 C and 80 C K values from
    Crowther 1979, which would give a **second temperature series** of the preheat effect; (ii) per-condition PDI
    and electrophoresis for the four autoclave states; (iii) the injected masses; (iv) an aqueous replication of
    the untreated-versus-autoclaved comparison, which is the single measurement that would make this paper's
    preheat factor shippable.
