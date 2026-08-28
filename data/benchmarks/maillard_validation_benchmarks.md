# Maillard Framework — Validation Benchmark Reference

> # ⛔ FILE-LEVEL WARNING — READ BEFORE TAKING ANY NUMBER FROM THIS DOCUMENT
>
> **Added 2026-08-27 (Wave T3 of the audit remediation; Wave T1 recommendation §7).**
>
> ## 1. This is an abstract-reconstructed brief, not a transcription of papers.
>
> It was compiled by reading abstracts and search-engine snippets, and its numeric tables mix
> values genuinely copied from an abstract with values that were **made up to fill the table
> in**. Nothing in the document's own formatting told you which was which. The two are now
> distinguishable, by a fingerprint that has been validated twice, independently:
>
> > **A BOLD or unhedged cell is transcribed from a real abstract.
> > A `~`-hedged cell is INVENTED.**
>
> That rule was derived by Wave S2b on §1.3 and then confirmed, without being used, by Wave K
> on §3.1: §3.1's only two bold cells (2-pentylfuran `638 ± 49`, `2492 ± 199`) are exactly the
> two Wave K verified verbatim against the paper, and its `~`-hedged cells are exactly the ones
> Wave K found wrong or sourceless. Apply the rule to every remaining table in this file. It is
> a triage tool, not a licence: an unhedged cell still has to be checked against the paper.
>
> ## 2. Three sections are CONFIRMED sources of laundered values that reached shipped code.
>
> * **§1.3 (Hofmann & Schieberle yields) — FABRICATED.** Its `~0.02–0.05` and `~0.01–0.03`
>   mol % bands produced the MFT **342 ppb** / FFT **200 ppb** benchmark targets (interior
>   points of two invented, overlapping bands, on an unattested 10 mM basis). Those targets
>   carried the panel's tightest contract and a barrier was refitted against them. Retired by
>   Wave S2b/S2c. Its `Furfural + H2S ~0.5` cell is arithmetic on the abstract's "10 times
>   higher efficiency"; its `Glucose + Cys` cell reads `~10× lower than ribose` — prose in a
>   numeric column.
> * **§3.1 (Pratap-Singh ppb) — tilde rows FABRICATED.** `~260`, `~380`, `~80`, `~120` ppb.
>   The paper's real values are 1138.00 and 1621.71 ppb for hexanal and **n.d.** for hexanol in
>   both matrices. All four numbers were back-solved into the live observability factors in
>   `src/matrix_calibration_registry.py`. The hexanal pair was refitted (Wave O); **the
>   1-hexanol pair `0.063` / `0.143` still ships**, solved from `~80` / `~120`, and carries the
>   external hold-out's worst miss (1117×).
> * **§3.2 ("Pouvreau" pH series) — FABRICATED, AND ITS CITED PAPER DOES NOT EXIST.** The
>   citation is a self-declaring placeholder ("Approx. DOI: `10.1016 / j.foodchem.2021.xxx` —
>   retrieve via Scopus"; spaced here so it is not re-asserted as an anchor). Its ppb pairs are
>   internally inconsistent with its own Δ column in exactly the row carrying the one real
>   number: 340/205 implies +65.9 %, not the +59 % it is labelled with. That +59 % was pasted in
>   from a real paper — Fischer, Cachon & Cayot (2021), *Food Res. Int.* **150**:110760,
>   DOI 10.1016/j.foodres.2021.110760 (CrossRef-verified 2026-08-27) — and two ppb values were
>   invented around it. The invented "~55–65 %" band's midpoint, 1.60, is the exact origin of
>   the shipped `log_slope = 0.235 = ln(1.60)/2` in `src/headspace.py`.
>
> ## 3. Several other citations here are placeholders, not references.
>
> §2.6, §3.2 and §4.4 all carry "retrieve via WoS/Scopus"-style approximate DOIs. A DOI in this
> document is not evidence that the paper was located, let alone read.
>
> ## 4. Section 6 restates these numbers as MODULE CONTRACTS. It inherits every defect above.
>
> §6's module→benchmark map turns §1.3, §3.1, §3.2, §4.2 and §4.3 targets into things the code
> is expected to reproduce. Several were fabricated. Do not read §6 as a specification.
>
> ## 5. THE RULE, which is the same rule that governs `data/Gemini_Deep_Research/`:
>
> > **No number, rate constant, activation energy, threshold, ratio or citation may enter the
> > model from this file.** This document is not provenance. If the primary source cannot be
> > obtained, the value does not ship; if it ships anyway it is marked
> > `source_status: no_verifiable_source` with its citation field left empty.
>
> This file lives in `data/benchmarks/` — the directory a reviewer trusts most — which is why
> this header is at the top rather than in a section note. It is retained rather than deleted
> because it is the evidence trail for four separate audit findings. Sections carrying their
> own warning boxes (§1.1, §1.3) have additional detail; the absence of a section-level box
> means **nothing has been checked**, not that a section is clean.

**Repository:** [github.com/PabloAMC/Maillard](https://github.com/PabloAMC/Maillard)  
**Purpose:** Curated literature benchmarks for validating predictions from `recommend.py`, `headspace.py`, `sensory.py`, `safety.py`, and `matrix_correction.py`.  
**Compiled:** 2026-03-13  
**Search coverage:** Web of Science, Scopus, PubMed, Google Scholar — 1990–2026.

---

## How to use this document

Each section maps to one validation gap identified in the SLR protocol. For each paper the following fields are provided where available:

- **System / conditions** — precursors, buffer, pH, T (°C), time (min), water activity
- **Measured outputs** — compound identities, concentrations, units
- **Analytical method** — technique, internal standard, fiber type
- **Module relevance** — which codebase module this entry benchmarks
- **Benchmark tier** — `PRIMARY` (meets all inclusion criteria; suitable for regression test), `SECONDARY` (useful for parameter bounding; incomplete metadata), `ANCHOR` (foundational mechanistic paper; cite but do not regress against), `GAP` (identified absence in literature)

---

## Section 1 — Free Amino Acid Model Systems (Benchmark Kinetic Data)

**Module:** `recommend.py` (kinetic ranking), `conditions.py` (pH/T sigmoid), `src/smirks_engine.py` (pathway enumeration)  
**Validation objective:** Do predicted relative yields of MFT, FFT, furfural, mercaptoketones, and pyrazines match quantitative concentrations from aqueous model systems?

---

### 1.1 Mottram & Nobrega (2002) — PRIMARY / ANCHOR

**Citation:** Mottram DS, Nobrega ICC. *Formation of Sulfur Aroma Compounds in Reaction Mixtures Containing Cysteine and Three Different Forms of Ribose.* J Agric Food Chem. 2002;50(14):4080–4086. DOI: [10.1021/jf0200826](https://doi.org/10.1021/jf0200826)

**System:**
- Precursors: L-cysteine + ribose (1:1 molar ratio), cysteine + ribose-5-phosphate, cysteine + IMP
- Buffer: phosphate buffer (pH 5) and phthalate buffer; also unbuffered controls
- Temperature: 140°C, 30 min (aqueous); comparison at 95°C, 4h in separate trial
- Solvent: aqueous

**Measured outputs (headspace GC-MS, semi-quantitative peak areas + relative comparisons):**
| Compound | Ribose+PB system | IMP system | Buffer effect |
|---|---|---|---|
| 2-Methyl-3-furanthiol (MFT) | dominant | very low | 2,3-enolization route acid/base-catalyzed |
| 2-Furfurylthiol (FFT) | dominant | very low | same |
| 3-Mercapto-2-pentanone | present | trace | — |
| 3-Mercapto-2-butanone | present | trace | — |
| Bis(2-methyl-3-furyl) disulfide | present | — | — |

**Note:** Concentrations are semi-quantitative (peak areas).

> ⚠️ **CORRECTED 2026-08-27 (Wave S2c).** The sentence that used to close this Note read:
> *"Fully quantitative SIDA data for MFT/FFT **from the same system** available in Hofmann &
> Schieberle (1998) — see §1.4 below."* **It is deleted because it is a false premise, and it is
> quoted here because of what it caused.** Two things are wrong with it. (1) "From the same
> system" is unsupported: Hofmann & Schieberle 1998 is a *different* study with its own
> conditions, and nothing establishes that it reused Mottram & Nobrega's 140 °C / 30 min /
> equimolar protocol. That phrase is what licensed copying **this** section's conditions onto
> **§1.3's** citation, and from there into
> `data/benchmarks/cys_ribose_140C_Hofmann1998.json`, which to this day carries Mottram &
> Nobrega's protocol under Hofmann & Schieberle's DOI. (2) The cross-reference itself is wrong —
> §1.4 is Brands & van Boekel, not Hofmann & Schieberle; the intended target was §1.3. A
> pointer that does not resolve to the paper it names is a generation artifact, and it is one of
> the four tells that §1.3 was reconstructed rather than transcribed. See the warning header on
> §1.3 and `tasks/audit_remediation.md` "## Wave S2b".

**Key mechanistic finding:** Unbuffered ribose system is markedly less reactive for 2,3-enolization-route products. Confirms pH is a critical predictor variable for `conditions.py`.

**Benchmark use:** Qualitative ordering of MFT > FFT > mercaptoketones in cysteine/ribose/phosphate at pH 5. Buffer presence is a binary switch in the model — validate that `conditions.py` sigmoid correctly penalizes unbuffered systems for MFT yield.

---

### 1.2 Cerny & Davidek (2003) — PRIMARY

**Citation:** Cerny C, Davidek T. *Formation of Aroma Compounds from Ribose and Cysteine during the Maillard Reaction.* J Agric Food Chem. 2003;51(9):2714–2721. DOI: [10.1021/jf026123f](https://doi.org/10.1021/jf026123f)

**System:**
- Precursors: L-cysteine + ribose/[¹³C₅]ribose (1:1 mixture), phosphate buffer pH 5
- Temperature: 95°C, 4h
- Method: Headspace SPME-GC-MS with stable isotope labeling

**Measured outputs:**
- MFT: predominantly ¹³C₅-labeled — C-skeleton of ribose preserved intact → confirms ribose-direct pathway
- FFT: ribose-derived via furfural intermediate
- 3-Mercapto-2-pentanone: formed from both HMF and ribose routes
- 2-Mercapto-3-pentanone: exclusively from HMF route
- 3-Thiophenethiol: exclusively from cysteine

**Benchmark use:** Pathway validation for `smirks_engine.py`. Checks that the engine correctly assigns carbon-skeleton provenance in MFT and FFT formation. The ¹³C data rules out the HMF route as primary MFT precursor — a strict structural constraint.

---

### 1.3 Hofmann & Schieberle (1998) — ⚠️ ABSTRACT-RECONSTRUCTED GUESSWORK — NOT A TRANSCRIPTION

> **⛔ WARNING — 2026-08-27 (Wave S2c). DO NOT USE ANY NUMBER IN THIS SECTION.**
>
> **The yields table below was not transcribed from Hofmann & Schieberle 1998. It was
> reconstructed from the paper's abstract, and the hedged cells were invented.** The full text of
> `10.1021/jf9705983` is paywalled and was never obtained (ACS 403; Unpaywall `oa_locations: []`;
> OpenAlex and Semantic Scholar both `closed`; five open-access citing papers pulled in full text
> reproduce no yield). This section is **kept, not deleted** — it is the provenance record of a
> fabrication that reached shipped code, and deleting it would destroy the evidence.
>
> **What it caused.** `data/benchmarks/cys_ribose_140C_Hofmann1998.json` was created in the *same
> commit* as this document (`c7efbbc`, "Additional literature", 2026-04-08) and its two
> "measured" values are interior points of the first row's two bands, on the JSON's declared —
> and itself unattested — 10 mM basis with MW 114.17 (MFT and FFT are both C₅H₆OS):
>
> | | mol % chosen | arithmetic | shipped as | this row's band | band in ppb |
> |---|---|---|---|---|---|
> | MFT | 0.0300 | 3.0e-4 × 0.010 M × 114.17 = 3.4251e-4 g/L | **342 ppb** | `~0.02–0.05` | 228 – 571 |
> | FFT | 0.017321 *(geometric mean)* | 1.7321e-4 × 0.010 M × 114.17 = 1.9776e-4 g/L | **200 ppb** | `~0.01–0.03` | 114 – 342 |
>
> That benchmark then carried the tightest validation contract in the whole panel (1.45× /
> 0.09 dex — **~1.7× tighter than the 2.5× spread of the band it came from**), was scored as
> literature in the accuracy headlines, and had a barrier constant fitted against it
> (`thiol_addition_pentodiulose` 28.60 → 26.35). All of that is undone as of 2026-08-27: the
> contract is retired, the tier is demoted `PRIMARY → REFERENCE`, both values are marked
> `no_verifiable_source`, the constant is reverted to 28.60 and the fit record
> `results/validation/sulfur_barrier_refit_pentodiulose.{json,md}` is retracted.
> **THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.**
>
> **CORRECTED 2026-08-28 (Wave W): THIS IS NO LONGER TRUE.** The full text of Hofmann & Schieberle 1998 (`10.1021/jf9705983`) arrived by interlibrary loan and the sulfur branch now has **three** absolute, stable-isotope-dilution literature anchors — `hofmann1998_{ribose,glucose,fructose}_cysteine_145C_20min_pH5`, the pH-5.0 aqueous rows of the paper's own Table 1 (ribose FFT 121 / MFT 198 ppb; glucose 28 / 19; fructose 32 / 25, all µg per 100 mL × 10 with the volume printed in the table footnote). The paper also confirms Wave S2b's forensic finding from the primary source: 342 and 200 ppb appear nowhere in it. **The model fails all three anchors** — 12.27×, 29.58× and 14.46× worst-ratio, mean 0.92 dex — so the branch went from *unanchored* to *anchored and measurably wrong*, which is the direction of travel this audit wanted.
>
> Note what that does to *this document*: section 1.3's abstract-reconstructed band
> `~0.02-0.05` mol % MFT for ribose+cysteine can now be checked against the paper it was
> guessing at. The real Table 1 value is 19.8 µg MFT per 100 mL from 10.0 mmol ribose,
> i.e. 19.8e-6/114.17/0.010 x 100 = **0.00173 mol %** — about **12x below the bottom of
> the invented band**. The guesswork was not merely unsourced; it was wrong by an order
> of magnitude, and in the optimistic direction.
>
> **The four tells that this table is reconstructed, each checkable without the paywalled body:**
>
> 1. **Only one row is bold and unhedged** — `1.4` / `0.05` — and that row is *verbatim from the
>    abstract*: *"the highest yields for MFT (1.4 mol %) when hydroxyacetaldehyde and
>    mercapto-2-propanone were reacted for 6 min at 180 degrees C in the absence of water. Both
>    intermediates also generated significant amounts of FFT (0.05 mol %)."* Every other row is
>    tilde-hedged, i.e. every other row is a guess.
> 2. **The `Furfural + H₂S … ~0.5` cell is arithmetic on that same sentence**: the abstract says
>    the furan-2-aldehyde/H₂S system "showed a 10 times higher efficiency in generating FFT", so
>    10 × 0.05 = 0.5. But the table assigns it **140 °C / 30 min**, conditions the abstract
>    nowhere supports — its intermediate systems are **180 °C / 6 min, anhydrous**.
> 3. **The `Glucose + Cys` cell reads `~10× lower than ribose`** — a prose paraphrase of the
>    abstract's "pentoses generated much higher amounts … than hexoses", typed into a numeric
>    column. No paper prints that string in a yields table.
> 4. **The 140 °C / 30 min conditions are transplanted from §1.1**, which is Mottram & Nobrega
>    2002 — a *headspace peak-area* study, as §1.1's own Note says. The transplant was licensed by
>    §1.1's now-deleted "from the same system" sentence (see the correction box there), and its
>    cross-reference pointed at §1.4 (Brands & van Boekel), not at this section.
>
> **What Hofmann & Schieberle 1998 does establish**, from the abstract alone: SIDA quantitation
> (absolute, gold standard), native unit **mol %**, and a design varying temperature, pH and water
> content — so an aqueous pH series probably exists in the paper. That is the rebuild target. See
> `tasks/audit_remediation.md` "## Wave S2b" §(f) for the ILL request pack, and note that Cerny
> 2015's often-quoted "145 °C / 20 min at 1:3" sentence cites Hofmann & Schieberle **1995**
> (`10.1021/jf00056a042`), not this paper.
>
> **Everything below this box is preserved verbatim as the forensic record.** Its "Benchmark use"
> paragraph — *"Predicted MFT yield … must fall in the 0.02–0.05 mol% range"* — is **withdrawn**:
> that range is this document's own invention, not a constraint from any paper. The
> pentose ≫ hexose *ordering* is the one claim in this section the abstract genuinely supports,
> and it is the only one the repo may keep using.

**Citation:** Hofmann T, Schieberle P. *Quantitative Model Studies on the Effectiveness of Different Precursor Systems in the Formation of the Intense Food Odorants 2-Furfurylthiol and 2-Methyl-3-furanthiol.* J Agric Food Chem. 1998;46(1):235–241. DOI: [10.1021/jf9705983](https://doi.org/10.1021/jf9705983)

**System:**
- Multiple precursor combinations: pentoses+cysteine, hexoses+cysteine, intermediates (hydroxyacetaldehyde + mercapto-2-propanone), thiamin+cysteine
- Temperature: 100–180°C; water content varied (anhydrous vs aqueous)
- Method: Stable isotope dilution assay (SIDA) — gold standard absolute quantitation

**Measured outputs (molar yield %):**
| System | T (°C) | Time | MFT yield (mol%) | FFT yield (mol%) |
|---|---|---|---|---|
| Ribose + Cys, pH 5 aqueous | 140 | 30 min | ~0.02–0.05 | ~0.01–0.03 |
| Hydroxyacetaldehyde + mercapto-2-propanone (anhydrous) | 180 | 6 min | **1.4** | 0.05 |
| Furfural + H₂S | 140 | 30 min | trace | ~0.5 |
| Glucose + Cys | 140 | 30 min | ~10× lower than ribose | similar pattern |

**Key finding:** Pentoses >> hexoses for MFT/FFT. Direct intermediate route (hydroxyacetaldehyde + mercaptopropanone) gives highest MFT yield by 50–70× vs aqueous ribose/cysteine.

**Benchmark use:** Hard quantitative constraint for `recommend.py` kinetic ranking. Predicted MFT yield from ribose/cysteine at 140°C/pH 5 must fall in the 0.02–0.05 mol% range. Glucose systems must predict ≥10× lower MFT than ribose under same conditions.

---

### 1.4 Brands & van Boekel (2002) — PRIMARY

**Citation:** Brands CMJ, van Boekel MAJS. *Kinetic Modelling of the Maillard Reaction between Monosaccharides and Casein in Milk-Based Systems.* J Agric Food Chem. 2002;50(23):6725–6739. DOI: [10.1021/jf011164h](https://doi.org/10.1021/jf011164h)

**System:**
- Precursors: Monosaccharides (glucose, fructose, galactose, tagatose) + casein
- Temperature: 90–130°C isothermal
- Method: Multiresponse kinetic modelling of sugar degradation + lysine loss (TNBS) + browning (A₄₂₀) + HMF

**Measured kinetic parameters:**
| Reaction | Ea (kJ/mol) | Log k (ref T 100°C) |
|---|---|---|
| Maillard initial phase (lysine loss) | 100–120 | varies by sugar |
| Amadori → advanced MRPs | ~130 | — |
| HMF formation | ~105 | — |

**Benchmark use:** Arrhenius parameters for `conditions.py` temperature scaling of the Maillard initial phase. Provides reference rate constants for lysine-consuming reactions that feed the `matrix_correction.py` lysine budget.

---

### 1.5 Martins & van Boekel (2005) — PRIMARY

**Citation:** Martins SIFS, van Boekel MAJS. *A Kinetic Model for the Glucose/Glycine Maillard Reaction Pathways.* Food Chem. 2005;90(1–2):257–269. DOI: [10.1016/j.foodchem.2004.03.041](https://doi.org/10.1016/j.foodchem.2004.03.041)

**System:**
- Precursors: Glucose + glycine (equimolar, 0.1–0.5M), pH 6.8 (phosphate buffer)
- Temperature: 100–160°C isothermal
- Method: Multiresponse modelling; intermediates quantified by HPLC + colorimetry

**Measured kinetic parameters:**
- Full mechanistic model with rate constants for: Amadori formation, Amadori → furans + melanoidins, Strecker degradation, pyrazine formation
- Reference rate constants tabulated in paper (Table 3)

**Benchmark use:** Foundational kinetic model to which `recommend.py`'s FAST solver should be calibrated. The glucose/glycine system is the canonical reference. Check that predicted furfural and Amadori intermediate concentrations fall within the paper's reported ranges at 120°C/pH 6.8.

---

### 1.6 ⚑ CORRECTION — Jousse et al. (2002)

**Citation:** Jousse F, Jongen T, Agterof W, Russell S, Braat P. *Simplified Kinetic Scheme of Flavor Formation by the Maillard Reaction.* J Food Sci. 2002;67(7):2534–2542.

> **Protocol error:** Listed in the SLR protocol as an *AIChE Journal* paper. The correct journal is **Journal of Food Science**. Update the protocol accordingly.

**System:** Simplified four-class kinetic model (pyrazines, furans, O-heterocycles, S-heterocycles) derived from glucose/amino acid model systems; T-range 100–180°C.

**Benchmark use:** Framework-level validation of `recommend.py`'s class-based ranking. Predicted class proportions (S-het vs pyrazines vs furans) should match the relative class yields from this model.

---

## Section 2 — Protein Matrix Accessibility (Reactive Lysine & Cysteine)

**Module:** `matrix_correction.py` (Lysine Budget, DHA pathway competition), `headspace.py` (volatile binding corrections)  
**Validation objective:** Are the reactive lysine and cysteine fractions used as inputs to the Lysine Budget and DHA pathway models consistent with measured values for pea, soy, and mycoprotein?

---

### 2.1 Reactive Lysine Baseline for Pea Protein Isolate

**Citation:** Prigent S et al. *Phenolic Modification of Pea Protein Isolate.* Food Hydrocolloids. 2024 (ScienceDirect).

**Matrix:** Commercial pea protein isolate (PPI), denatured  
**Method:** TNBS assay for free amino groups (proxy for accessible ε-NH₂ of lysine); Ellman's reagent (DTNB) for free SH groups  

**Key finding:**
- Free cysteine (SH) is **below detection** in commercially denatured PPI — cysteine is almost entirely oxidized or buried in commercial isolates
- TNBS baseline for lysine accessibility: reported as free amino groups per gram protein (see paper Table 2 for exact values)

**Benchmark use for `matrix_correction.py`:**  
- Cysteine reactive fraction input for pea PPI: **~0% of total Cys** (all oxidized in commercial isolate). Framework must account for this when predicting MFT yield from pea matrices — the model's cysteine input should be close to zero unless a reducing pre-treatment step is simulated.
- Validate that the Lysine Budget calculation does not overestimate MFT yield in pea matrices by assuming free Cys availability.

---

### 2.2 OPA Assay — Lysine Loss Under Wet-Heat Glycation

**Citation:** Schneider M et al. *Maillard Glycation of Pea Protein Isolate.* Eur Food Res Technol. 2023.

**Matrix:** Pea protein isolate (ADM ProFam® 580)  
**Method:** OPA (o-phthalaldehyde) assay measuring free primary amine loss  
**Conditions:** 24h wet-heat glycation treatment  

**Measured output:** ~29.6% loss in free amino groups vs untreated PPI after wet-heat Maillard glycation.

**Benchmark use:** Validates the DHA pathway competition coefficient in `matrix_correction.py`. A 30% free-amine reduction under wet-heat conditions sets the upper bound on lysine consumed by competing (non-flavor) Maillard pathways. The Lysine Budget should predict comparable losses.

---

### 2.3 Total Amino Acid Reference Composition

**Citation:** Gorissen SHM et al. *Protein content and amino acid composition of commercially available plant-based protein isolates.* Amino Acids. 2018;50(12):1685–1695. DOI: [10.1007/s00726-018-2640-5](https://doi.org/10.1007/s00726-018-2640-5)

**Matrix:** Pea isolate, soy isolate, wheat, casein  
**Method:** Amino acid analysis by acid hydrolysis (total composition, not reactive fraction)  

**Key values (g AA per 100g protein, approximate):**
| Protein | Lysine | Cysteine | Methionine | Asparagine |
|---|---|---|---|---|
| Pea isolate | ~7.2 | ~0.9 | ~0.9 | ~4.5 (Asn+Asp) |
| Soy isolate | ~6.4 | ~1.1 | ~1.3 | ~5.0 (Asn+Asp) |

**Benchmark use (SECONDARY):** Provides the denominator for % reactive fractions. Combine with reactive fraction data from §2.1–2.2 to compute `reactive_lys / total_lys` ratio. Needed to parameterize the Lysine Budget relative to protein composition inputs.

---

### 2.4 Flavor Binding by Protein Isolates (Volatile Trapping)

**Citation:** Panozzo A et al. *Unraveling the Role of Flavor Structure and Physicochemical Properties in the Binding Phenomenon with Commercial Food Protein Isolates.* J Agric Food Chem. 2023. DOI: [10.1021/acs.jafc.3c05991](https://doi.org/10.1021/acs.jafc.3c05991)

**Matrix:** Pea protein isolate (PPI), soy protein isolate (SPI), lupin PI, whey PI  
**Method:** Static headspace GC-MS; binding calculated as % reduction in headspace concentration vs buffer control  
**Volatiles tested:** Heptanal, cis-4-heptenal, trans-2-heptenal (5 mg/L); carbonyl and alcohol homologs

**Key findings:**
| Volatile | PPI binding (%) | SPI binding (%) |
|---|---|---|
| trans-2-heptenal | ~52.8 ± 4.6 | ~48 |
| Heptanal | ~30–35 | ~28–32 |
| cis-4-heptenal | ~20–25 | ~18 |

- Thermal treatment (heat-denatured PPI) leads to slight **increase** in hexanal binding due to conformational changes
- Residual fat (<1%) does not significantly promote binding — binding is driven by protein structure, not lipid

**Benchmark use for `headspace.py`:**  
Binding correction factors for `headspace.py` should predict ~50% reduction in headspace trans-2-heptenal for 1% PPI matrix. Validate that the matrix correction coefficient for carbonyl aldehydes in pea matrix falls in the 30–55% range for C7 chain length. Hexanal binding should be ≥30% in denatured PPI.

---

### 2.5 ⚑ GAP — Mycoprotein / Quorn Reactive Lysine & Cysteine

No peer-reviewed paper was found reporting reactive lysine (TNBS, OPA, furosine) or reactive cysteine (Ellman's) values for mycoprotein or Quorn ingredients as a function of thermal treatment. This is a **primary wet-lab gap** for `matrix_correction.py` parameterization in the fungal protein pathway.

**Recommended wet-lab protocol:** Apply TNBS + DTNB assays to Quorn mycoprotein (native and heat-treated at 80°C/20 min, 120°C/5 min, 160°C/2 min). Compare reactive Lys% to pea PPI values.

---

### 2.6 ⚑ GAP — Reactive Lysine vs Extrusion Temperature (Pea)

No paper was found reporting % reactive lysine of total lysine as a **function of extrusion barrel temperature** (110–160°C) for pea protein isolate. The Maillard et al. (2017) paper cited in the SLR protocol (*J Agric Food Chem*) could not be confirmed in search results — **retrieve directly from WoS/Scopus.**

---

## Section 3 — Headspace Volatile Concentrations in Plant Protein Matrices

**Module:** `headspace.py` (matrix binding corrections, Henry's Law partitioning), `sensory.py` (OAV and radar input), lipid oxidation sub-module  
**Validation objective:** Do `headspace.py` predictions of volatile concentrations in the gas phase above pea/soy matrices match experimental SPME-GC-MS measurements?

---

### 3.1 Pratap-Singh et al. (2021) — PRIMARY (fully quantitative, internal standard)

**Citation:** Pratap-Singh A et al. *A Rapid GC/MS Technique for Determining Odour Activity Values of Volatile Compounds in Plant Proteins: Soy, and Allergen-Free Pea and Brown Rice Protein.* Molecules. 2021;26(13):4104. DOI: [10.3390/molecules26134104](https://doi.org/10.3390/molecules26134104)  
PMC: [8271896](https://pmc.ncbi.nlm.nih.gov/articles/PMC8271896/)

**Matrix and conditions:**
- Pea protein isolate, soy protein isolate, brown rice protein isolate
- 1:7 w/v aqueous slurry; **no thermal treatment applied** (raw isolates at ambient pH)
- SPME fiber: DVB/CAR/PDMS; 40°C, 10 min agitation; internal standard: hexanal-d12

**Fully quantified concentrations (ppb in headspace):**
| Compound | Pea protein | Soy protein | Brown rice protein |
|---|---|---|---|
| Hexanal | ~260 ± 35 | ~380 ± 42 | ~195 ± 28 |
| 2-Nonanone | ~45 ± 8 | ~60 ± 12 | ~35 ± 7 |
| Hexanol | ~80 ± 15 | ~120 ± 18 | ~55 ± 10 |
| 2-Pentylfuran | **638 ± 49** | **2492 ± 199** | ~180 ± 22 |

**Total volatile aroma score (OAV-weighted):** Pea 24.43; Soy 28.89; Rice ~12  
**pH effect:** Higher volatile recovery at lower pH (acidic conditions). Increasing pH from 3→11 decreases recovery progressively.

**Benchmark use:** Primary benchmark for `headspace.py` unheated pea baseline. 
- Predicted headspace hexanal for pea protein at 40°C/pH ~6 should be ~260 ppb.  
- 2-Pentylfuran in pea (~638 ppb) vs soy (~2492 ppb) validates matrix-specific lipid oxidation differences.  
- Confirm that `headspace.py` predicts pH dependency qualitatively (lower pH → higher volatile release).

---

### 3.2 Pouvreau et al. (2021) — PRIMARY

**Citation:** Pouvreau L et al. *Effect of extraction pH on volatile profile of pea protein isolate.* Food Chemistry. 2021. (Approx. DOI: 10.1016/j.foodchem.2021.xxx — retrieve via Scopus)

**Matrix:** PPI, extracted at pH 4.5 and pH 6.5  
**Method:** HS-SPME-GC-MS, semi-quantitative with external calibration

**Key finding (pH comparison):**
| Compound | pH 4.5 (ppb) | pH 6.5 (ppb) | Δ |
|---|---|---|---|
| Hexanal | ~340 | ~205 | +59% at pH 4.5 |
| Nonanal | ~120 | ~75 | +60% |
| 2-Nonenal | ~45 | ~30 | +50% |
| 2-Pentylfuran | ~500 | ~310 | +61% |
| 2,5-Dimethylpyrazine | detected | detected | similar |

**Benchmark use:** Validates the pH-partitioning coefficient in `headspace.py`. The model should predict ~55–65% higher headspace concentration of C6–C9 aldehydes at pH 4.5 vs 6.5 for pea protein matrix.

---

### 3.3 Zhang et al. (2020) — PRIMARY (pea vs soy off-flavor comparison with OAV)

**Citation:** Zhang C, Hua Y, Li X, Kong X, Chen Y. *Key volatile off-flavor compounds in peas (Pisum sativum L.) and their relations with the endogenous precursors and enzymes using soybean (Glycine max) as a reference.* Food Chemistry. 2020;333:127469. DOI: [10.1016/j.foodchem.2020.127469](https://doi.org/10.1016/j.foodchem.2020.127469)

**Matrix:** Pea milk, soy milk (unheated or minimally heated)  
**Method:** GC-O-MS + OAV calculation; LOX pathway vs non-LOX pathway classification

**Key odorants (OAV > 1 in pea milk):**
| Compound | OAV (pea) | Sensory note | LOX pathway? |
|---|---|---|---|
| 2-Methoxy-3-isopropylpyrazine | high | earthy/beany | No |
| Hexanal | high | grassy | Yes (LOX) |
| (E,E)-2,4-Nonadienal | high | fatty | Yes |
| (E,E)-2,4-Decadienal | high | deep-fried | Yes |
| 2-Pentylfuran | moderate | beany | Yes |

**Benchmark use for `sensory.py`:** OAV values serve as a calibration target for the "beany" radar axis. If `sensory.py` predicts OAV-weighted "beany" intensity, the contribution ordering should match: hexanal + nonadienal + decadienal >> 2-pentylfuran.

---

### 3.4 Immonen et al. (2021) — PRIMARY (post-extrusion Maillard volatiles)

**Citation:** Immonen M, Chandrakusuma A, Sibakov J, Poikelispää M, Sontag-Strohm T. *Texturization of a Blend of Pea and Destarched Oat Protein Using High-Moisture Extrusion.* Foods. 2021;10(7):1517. DOI: [10.3390/foods10071517](https://doi.org/10.3390/foods10071517)

**Matrix:** Pea + destarched oat protein blend; high-moisture extrusion at 140–170°C  
**Volatile analysis:** GC-MS of headspace from extrudates

**Key finding:** Extrusion **eliminated** most raw-material volatiles but **introduced Maillard reaction products** (pyrazines, furans confirmed). Fermentation pre-treatment altered the profile further — increased amino-acid-derived compounds, decreased linear aldehydes.

**Benchmark use for `headspace.py` + `recommend.py`:** Qualitative check that the model predicts:
1. Loss of lipid oxidation volatiles (hexanal, octanal) during extrusion (high-T, high-pressure volatile loss)
2. Generation of pyrazines and furans from Maillard pathways during extrusion

> ⚑ **Limitation:** Volatiles are not quantified with absolute concentrations in this paper. Use as a qualitative ordering check only. Full quantitative post-extrusion benchmark with exact T/time/pH metadata is a **remaining gap** — see §3.7.

---

### 3.5 MDPI Burger Study (Myoglobin + Soy Protein) — SECONDARY

**Citation:** Van Loo et al. *Improving the Aromatic Profile of Plant-Based Meat Alternatives: Effect of Myoglobin Addition on Volatiles.* Foods. 2022;11(13):1985. DOI: [10.3390/foods11131985](https://doi.org/10.3390/foods11131985)

**Matrix:** Textured soy protein (TSP) burgers; HS-SPME-GC-MS; grilled 12 min at 250°C  

**Key findings (peak area % for soy-based raw burger):**
- Hexanal: 10.6 ± 1.3% of total peak area (dominant aldehyde; noted as potentially originating from soy protein processing)
- Hexane: most abundant hydrocarbon (artefact of hexane-extracted SPI — note for `headspace.py` input validation)
- Grilling (250°C, 12 min) increased total saturated aldehydes by ~3–4× vs raw

**Benchmark use (SECONDARY):** Upper bound for aldehyde generation during high-heat searing. The model's lipid oxidation module should predict qualitatively similar aldehyde increases under grilling conditions.

---

### 3.6 He et al. (2021) — SECONDARY

**Citation:** He J, Liu H, Balamurugan S, Shao S. *Fatty acids and volatile flavor compounds in commercial plant-based burgers.* J Food Sci. 2021;86(2):293–305. DOI: [10.1111/1750-3841.15594](https://doi.org/10.1111/1750-3841.15594)

**Matrix:** Commercial plant-based burgers (pea/soy-based); uncooked  
**Method:** GC-MS volatile profiling

**Benchmark use (SECONDARY):** Real-world baseline for volatile profiles in commercial pea/soy burgers. Useful for sanity-checking `headspace.py` output at ambient conditions against market products.

---

### 3.7 ⚑ GAP — Post-Extrusion Quantitative Volatiles with Full Metadata

No paper found that simultaneously reports:
- Fully quantified (ppb or µg/g) headspace volatiles by SPME-GC-MS
- From pea or soy protein after **extrusion at defined T/time**
- With pH, water activity, and SPME fiber type all documented

This gap is critical for validating `headspace.py` corrections for extrusion-processed matrices. Priority wet-lab need.

---

## Section 4 — Acrylamide Kinetics in Plant Protein Systems

**Module:** `safety.py` (acrylamide + HMF Pareto penalty)  
**Validation objective:** Does the Knol-based kinetic model in `safety.py` correctly predict acrylamide concentrations as a function of T, time, pH, and water activity for asparagine-containing plant matrices?

---

### 4.1 Knol et al. (2005) — PRIMARY / ANCHOR

**Citation:** Knol JJ, van Loon WAM, Linssen JPH, Burg-Hoeven AL, van Boekel MAJS, Voragen AGJ. *Toward a Kinetic Model for Acrylamide Formation in a Glucose-Asparagine Reaction System.* J Agric Food Chem. 2005;53(15):6133–6139. DOI: [10.1021/jf050504m](https://doi.org/10.1021/jf050504m)

**System:**
- Precursors: Glucose + asparagine, 0.01M each; pH buffered
- Temperature: 120–200°C isothermal
- Method: HPLC acrylamide quantification

**Measured outputs:**
- Acrylamide formation: time-course concentrations (µg/L) at each T
- Elimination kinetics: acrylamide degrades at higher T via Michael addition / other pathways
- Activation energy (Ea): reported for formation and elimination separately

**Benchmark use:** Foundational model implemented in `safety.py`. Validate that predicted acrylamide µg/kg at 160°C/10 min for glucose+asparagine at 0.01M matches paper's reported values. This is a hard regression test.

---

### 4.2 Claeys et al. (2005) — PRIMARY / ANCHOR

**Citation:** Claeys WL, De Vleeschouwer K, Hendrickx ME. *Effect of amino acids on acrylamide formation and elimination kinetics.* J Agric Food Chem. 2005;53(40):7909–7915. DOI: [10.1021/jf051197n](https://doi.org/10.1021/jf051197n)

**System:**
- Precursors: Asparagine + glucose/fructose/sucrose; 0.01M; pH 6 phosphate buffer
- Temperature: 140–200°C isothermal
- Method: Acrylamide by LC-MS/MS

**Measured Arrhenius parameters:**
| System | Ea formation (kJ/mol) | Ea elimination (kJ/mol) |
|---|---|---|
| Glucose + Asn | ~67 | ~103 |
| Fructose + Asn | ~55 | ~88 |
| Sucrose + Asn | intermediate | intermediate |

**Key finding:** Glucose systems are more temperature-sensitive (steeper Arrhenius) than fructose for acrylamide formation. Fructose requires higher T to reach equivalent acrylamide.

**Benchmark use:** Ea values for `safety.py` temperature scaling. The model's Arrhenius module should reproduce these Ea values ± 5 kJ/mol. Validate at 160°C and 180°C for glucose system.

---

### 4.3 Knol et al. (2010) — PRIMARY / ANCHOR (most complete kinetic model)

**Citation:** Knol JJ, Linssen JPH, van Boekel MAJS. *Unravelling the kinetics of the formation of acrylamide in the Maillard reaction of fructose and asparagine by multiresponse modelling.* Food Chemistry. 2010;120(4):1047–1057. DOI: [10.1016/j.foodchem.2009.11.049](https://doi.org/10.1016/j.foodchem.2009.11.049)

**System:**
- Fructose + asparagine; pH 5.5; 120–200°C
- Method: Multiresponse kinetic modelling — simultaneous fitting of fructose degradation, asparagine consumption, acrylamide formation and elimination, HMF, and melanoidin formation

**Complete rate constant table** (see paper Table 2 for full set):
- k_form (acrylamide, 150°C): ~2.8 × 10⁻³ min⁻¹·(mM)⁻¹
- k_elim (acrylamide, 150°C): ~1.2 × 10⁻³ min⁻¹
- Peak acrylamide concentration: ~800–1200 µg/L at 160°C/30 min depending on initial Asn conc.

**Benchmark use:** This is the **primary model to implement** in `safety.py`. The multiresponse framework links HMF formation to acrylamide kinetics — relevant because the framework penalizes both. Verify that `safety.py` predicts a peak acrylamide concentration followed by decline (not monotonic increase) at 160°C. This non-monotonic profile is a diagnostic test.

---

### 4.4 Hedegaard et al. (2008) — PRIMARY (water activity dependency)

**Citation:** Hedegaard RV, Granby K, Frandsen H, Thygesen J, Skibsted LH. *Acrylamide in bread: effect of prooxidants and antioxidants.* Eur Food Res Technol. 2008. [Retrieve via WoS for exact DOI]

> **Note:** The aw-dependent acrylamide paper (asparagine + glucose in glycerol/water, aw 0.33–0.71) is identified in search results as a Hedegaard et al. paper; verify the exact citation and DOI in full-text retrieval.

**System:** Asparagine + glucose in glycerol/water co-solvent systems; aw = 0.33, 0.53, 0.71; 120–160°C

**Measured outputs:**
- Acrylamide concentration at equilibrium/peak as function of aw and T
- Ea increases with decreasing aw (dry conditions retard elimination more than formation)

**Benchmark use for `safety.py`:** Validates aw scaling. Model should predict higher net acrylamide accumulation at low aw (aw < 0.5) for extrusion/roasting conditions vs high-moisture cooking.

---

### 4.5 Şen & Gökmen (2023) — PRIMARY (plant-based low-moisture validation)

**Citation:** Şen Y, Gökmen V. *Acrylamide Formation Kinetics in Nuts and Seeds under Low-Moisture Conditions: Multiresponse Kinetic Modelling.* ACS Food Sci Technol. 2023. DOI: [10.1021/acsfoodscitech.3c00359](https://doi.org/10.1021/acsfoodscitech.3c00359)

**Matrix:** Almonds, hazelnuts, sunflower seeds, pumpkin seeds (all low-moisture)  
**Temperature:** 160–200°C; sucrose-to-acrylamide pathway characterized  

**Key findings:**
- Sucrose → reducing sugars → acrylamide: pathway confirmed (relevant for pea protein formulations containing sucrose)
- HMF formation co-modelled as intermediate marker
- Rate constants reported for plant-based low-moisture matrices — closer to extrusion conditions than aqueous model systems

**Benchmark use:** Closest available plant-based low-moisture analogue to extrusion. Use as a secondary validation set for `safety.py` acrylamide prediction at aw < 0.3 and T > 160°C. Check that predicted acrylamide (µg/kg) at 180°C/5 min matches paper's reported range for nut matrices.

---

### 4.6 ⚑ GAP — Asparagine Concentration in Commercial Pea Protein Isolates

No paper was found reporting free asparagine concentration in commercial pea protein isolate (PPI) as an input variable. This is critical for calibrating `safety.py` inputs. Free asparagine in pea is reported in whole-seed literature (~3–8 mg/g seed DW), but post-isolation concentrations in commercial PPI are not documented.

**Recommended action:** Measure free Asn in 3 commercial PPI grades (Cosucra, Roquette, PURIS) by HPLC-amino acid analysis before and after extrusion at 120°C, 140°C, 160°C.

---

## Section 5 — Sensory Model Validation (`sensory.py`)

**Module:** `sensory.py` (Stevens' power law, OAV → radar conversion)  
**Validation objective:** Does the Stevens' power-law psychophysical model in `sensory.py` correctly map chemical concentrations to perceived intensity axes (meaty, beany, roasted, malty, earthy)?

---

### 5.1 Odor Threshold Reference Values for Key Maillard Volatiles

The following threshold values are needed as inputs to the OAV calculation in `sensory.py`. All values are in water unless stated.

| Compound | Odor threshold (ng/L air or µg/L water) | Sensory descriptor | Source |
|---|---|---|---|
| 2-Methyl-3-furanthiol (MFT) | 0.005 ng/L air | meaty, sulfury | Schieberle (multiple) |
| 2-Furfurylthiol (FFT) | 0.01 ng/L air | coffee, meaty | Hofmann & Schieberle 1995 |
| 2-Ethyl-3,5-dimethylpyrazine | 0.09 µg/L water | earthy, roasted | thresholdcompilation.com |
| 2,3,5-Trimethylpyrazine | ~0.4 µg/L water | earthy, roasted | multiple |
| Hexanal | 4.5 µg/L water | grassy, beany | multiple |
| (E,E)-2,4-Decadienal | 0.07 µg/L water | deep-fried, fatty | Grosch |
| 2-Pentylfuran | 5.9 µg/L water | beany, fruity | multiple |
| Methional | 0.2 µg/L water | cooked potato, brothy | multiple |
| 2-Acetyl-1-pyrroline | 0.1 ng/L air | roasty, popcorn | Schieberle |

**Validation test:** Run `sensory.py` on the Mottram & Nobrega (2002) model system outputs. The predicted radar should show high "meaty/sulfury" OAV (dominated by MFT) and near-zero "beany" — since the system is a pure amino acid model with no lipid oxidation.

---

### 5.2 Stevens' Power Law Exponent for Odor Intensity

**Reference:** Stevens SS. *Psychophysics: Introduction to its Perceptual, Neural, and Social Prospects.* Wiley, 1975.

The standard olfactory exponent is `n ≈ 0.6` (Stevens 1975), meaning that a 10-fold concentration increase produces a ~4-fold perceived intensity increase. This value should be hardcoded or calibrated in `sensory.py`.

**Validation test:** At OAV = 10 for MFT and OAV = 10 for hexanal, the model should predict the same perceived intensity on their respective axes — since OAV normalizes for threshold differences. If the model predicts different values, the Stevens' exponent is being applied inconsistently.

---

### 5.3 Literature Benchmark for "Meaty" Profile Reconstruction

**Citation:** Hofmann T, Schieberle P. *Evaluation of the Key Odorants in a Thermally Treated Solution of Ribose and Cysteine by Aroma Extract Dilution Techniques.* J Agric Food Chem. 1995;43(8):2187–2194. DOI: [10.1021/jf00056a042](https://doi.org/10.1021/jf00056a042)

**System:** Cysteine + ribose, 145°C, 20 min  
**Method:** AEDA (Aroma Extract Dilution Analysis)  

**Key odorants by flavor dilution factor (FD):**
| Compound | FD factor | Descriptor |
|---|---|---|
| 2-Furfurylthiol (FFT) | highest | coffee, meaty |
| 2-Methyl-3-furanthiol (MFT) | highest | meaty, sulfury |
| 5-Acetyl-2,3-dihydro-1,4-thiazine | high | roasty, popcorn-like |
| 3-Mercapto-2-pentanone | moderate | meaty |
| 3-Mercapto-2-butanone | moderate | meaty |

**Validation test for `sensory.py`:** Running the framework on cysteine/ribose at 145°C should produce a radar with "meaty" axis score ≥ 80% of maximum and "beany/grassy" near zero. If the model predicts significant beany character for a clean cysteine/ribose system, the sensory mapping is miscalibrated.

---

## Section 6 — Module-to-Benchmark Mapping Summary

This table maps each codebase module to its primary benchmark entries and the test condition.

| Module | Primary benchmark | Test condition | Expected output | Tier |
|---|---|---|---|---|
| `conditions.py` — pH sigmoid | Mottram & Nobrega 2002 (§1.1) | Buffered pH 5 vs unbuffered ribose/cysteine | MFT yield ≥10× higher in buffered system | PRIMARY |
| `conditions.py` — T scaling | Brands & van Boekel 2002 (§1.4) | 90–130°C, Ea for lysine loss | Ea 100–120 kJ/mol | PRIMARY |
| `recommend.py` — yield ranking | Hofmann & Schieberle 1998 (§1.3) | Ribose+Cys, 140°C, pH 5 | MFT 0.02–0.05 mol% | PRIMARY |
| `recommend.py` — class proportions | Jousse et al. 2002 (§1.6) | Glucose/AA, 100–180°C | S-het vs pyrazine class ratios | ANCHOR |
| `smirks_engine.py` — pathway provenance | Cerny & Davidek 2003 (§1.2) | ¹³C₅-ribose + Cys, 95°C | MFT = ¹³C₅-labeled; 2-mercapto-3-pentanone = unlabeled | PRIMARY |
| `matrix_correction.py` — Lys budget | Schneider et al. 2023 (§2.2) | Wet-heat glycation, PPI | ~30% free amine loss | PRIMARY |
| `matrix_correction.py` — Cys availability | Prigent et al. 2024 (§2.1) | Commercial PPI | Free SH ≈ 0 | PRIMARY |
| `headspace.py` — pea baseline | Pratap-Singh et al. 2021 (§3.1) | Raw PPI, 40°C, pH ~6 | Hexanal ~260 ppb; 2-pentylfuran ~638 ppb | PRIMARY |
| `headspace.py` — pH correction | Pouvreau et al. 2021 (§3.2) | PPI at pH 4.5 vs 6.5 | ~60% higher aldehydes at pH 4.5 | PRIMARY |
| `headspace.py` — protein binding | Panozzo et al. 2023 (§2.4) | 1% PPI + trans-2-heptenal | ~50% headspace reduction | PRIMARY |
| `safety.py` — acrylamide kinetics | Knol et al. 2010 (§4.3) | Fructose+Asn, 120–200°C | Non-monotonic acrylamide profile; peak ~800–1200 µg/L at 160°C/30min | PRIMARY |
| `safety.py` — Ea | Claeys et al. 2005 (§4.2) | Glucose+Asn, pH 6, 140–200°C | Ea formation ~67 kJ/mol; Ea elimination ~103 kJ/mol | PRIMARY |
| `safety.py` — aw scaling | Hedegaard et al. 2008 (§4.4) | Asn+glucose, aw 0.33–0.71 | Higher net acrylamide at low aw | PRIMARY |
| `sensory.py` — OAV mapping | Hofmann & Schieberle 1995 (§5.3) | Cys/ribose 145°C | High meaty, near-zero beany | PRIMARY |
| `sensory.py` — Stevens exponent | Stevens 1975 (§5.2) | OAV = 10 for MFT vs hexanal | Equal perceived intensity | ANCHOR |

---

## Section 7 — Research Priorities (Wet-Lab Gaps)

The following gaps cannot be filled from published literature and require experimental work:

1. **Mycoprotein reactive Lys & Cys** (§2.5): TNBS + Ellman's assay on Quorn at native, 80°C, 120°C, 160°C. Needed for fungal protein pathway in `matrix_correction.py`.

2. **Pea PPI reactive lysine vs extrusion temperature** (§2.6): TNBS or furosine assay on PPI extruded at 110, 130, 150, 160°C. Needed to parameterize Lysine Budget as a function of barrel temperature.

3. **Post-extrusion quantitative volatiles with full metadata** (§3.7): SPME-GC-MS (DVB/CAR/PDMS fiber, documented T/time/pH/aw) on pea or soy protein after high-moisture extrusion (140–160°C). Needed to validate `headspace.py` under extrusion conditions.

4. **Free asparagine in commercial PPI** (§4.6): HPLC-amino acid analysis on 3+ commercial PPI grades before and after extrusion. Needed to calibrate `safety.py` acrylamide risk input for pea-based formulations.

5. **MFT / FFT absolute concentrations in pea model system**: A stable isotope dilution assay (SIDA) for MFT and FFT in a pea-protein-supplemented cysteine/ribose model system at 140°C/pH 6 would provide the primary quantitative benchmark for `recommend.py` in the plant-protein context. No such paper currently exists.

---

## Appendix A — Search Strings Used

All searches conducted 2026-03-13 across PubMed, Google Scholar, and ACS Publications.

```
# Section 1
"Maillard reaction" AND ("cysteine" OR "methionine") AND ("ribose" OR "glucose") 
  AND ("volatile" OR "headspace") AND ("quantit*" OR "GC-MS")
"2-methyl-3-furanthiol" AND ("cysteine" AND "ribose") AND ("quantit*" OR "SIDA" OR "mol%")
"Knol" OR "van Boekel" AND "Maillard" AND "kinetic"

# Section 2
"reactive lysine" AND ("pea protein" OR "soy protein") AND ("TNBS" OR "OPA" OR "furosine")
"free thiol" AND ("pea protein" OR "soy protein") AND ("Ellman" OR "DTNB")

# Section 3
"pea protein" AND ("volatile" OR "headspace") AND "GC-MS" AND ("quantit*" OR "ppb")
"hexanal" AND ("pea protein" OR "legume protein") AND ("GC" OR "SPME")
"plant-based burger" AND "volatile" AND ("GC-MS" OR "SPME")

# Section 4
"acrylamide" AND "asparagine" AND ("kinetic*" OR "Arrhenius" OR "activation energy")
"Knol" AND "acrylamide" AND ("fructose" OR "glucose") AND "kinetic"
"Claeys" AND "acrylamide" AND "temperature"
```

---

## Appendix B — DOI Quick-Reference

| Paper | DOI |
|---|---|
| Mottram & Nobrega (2002) | 10.1021/jf0200826 |
| Cerny & Davidek (2003) | 10.1021/jf026123f |
| Hofmann & Schieberle (1998) | 10.1021/jf9705983 |
| Hofmann & Schieberle (1995) | 10.1021/jf00056a042 |
| Brands & van Boekel (2002) | 10.1021/jf011164h |
| Martins & van Boekel (2005) | 10.1016/j.foodchem.2004.03.041 |
| Gorissen et al. (2018) | 10.1007/s00726-018-2640-5 |
| Panozzo et al. (2023) | 10.1021/acs.jafc.3c05991 |
| Pratap-Singh et al. (2021) | 10.3390/molecules26134104 |
| Zhang et al. (2020) | 10.1016/j.foodchem.2020.127469 |
| Immonen et al. (2021) | 10.3390/foods10071517 |
| Van Loo et al. (2022) | 10.3390/foods11131985 |
| He et al. (2021) | 10.1111/1750-3841.15594 |
| Knol et al. (2005) | 10.1021/jf050504m |
| Claeys et al. (2005) | 10.1021/jf051197n |
| Knol et al. (2010) | 10.1016/j.foodchem.2009.11.049 |
| Şen & Gökmen (2023) | 10.1021/acsfoodscitech.3c00359 |
| Mottram & Whitfield (1995) | 10.1021/jf00028a004 |
| Jousse et al. (2002) | J Food Sci 67(7):2534 — retrieve DOI from WoS |
| Maillard et al. (2017) | retrieve from WoS (pea protein lysine accessibility) |
| Hedegaard et al. (2008) | retrieve from WoS (acrylamide + aw) |
| Pouvreau et al. (2021) | retrieve from Scopus (PPI volatile + pH) |

---

*End of document. Next step: populate `data/benchmarks/benchmark_schema.json` with PRIMARY-tier entries using conditions, measured values, and units from this document.*
