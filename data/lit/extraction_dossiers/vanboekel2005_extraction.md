# van Boekel — POSTER ABSTRACT, "Kinetic Modelling of the Maillard Reaction between Glucose and Glycine" (one page, p. 451 of *The Maillard Reaction in Foods and Medicine*; announces the study the trunk is fitted on, BEFORE the model existed)

### THE TRUNK'S OWN STUDY, ANNOUNCED AS WORK IN PROGRESS: this one-page abstract contains no rate constant, no barrier and no concentration — its last paragraph is "Research is underway to propose a kinetic model" — but it prints three conditions that DIFFER from the `MARTINS_M4` block the engine ships (pH **7** not 6.8; **40**-120 °C not 80-120 °C; and the Amadori compound named **fructosyl-lysine**, three times, in a glucose-and-glycine system whose Amadori compound is fructosyl-glycine).

**Source on disk:** `data/articles/vanboekel2005.pdf` (**one page**, 432 × 666 pt, 451 kB).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/vanboekel2005.txt`, 54 lines) **and then re-read from a 200 dpi render of the
page** (`scratchpad/vb2005-1.png`), because the text layer is an OCR of a scan with visible
corruption — "Arnadori" for "Amadori", "Chemise" for "Chemistry", "Nerh." for "Neth.",
"WageningenAgricultural", "The NetherlanA", "complicatesthe", "sofar". **Every quotation and every
number below was verified against the page image**, which is clean and fully legible. The page also
carries a second, unrelated abstract (van Chuyen et al. on dietary antioxidants and diabetes
mellitus); it is not extracted here beyond noting that it shares the page.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetic Modelling of the Maillard Reaction between Glucose and Glycine" |
| Author | **Martinus A. J. S. van Boekel** — sole author. "Wageningen Agricultural University, Department of Food Science, PO Box 8129, 6700 EV, Wageningen, The Netherlands" |
| Type | **Poster abstract**, one page, under the running head "Poster Abstracts", page number **451** |
| Book (from the PDF's own metadata) | Subject field: **"The Maillard Reaction in Foods and Medicine (2005) 451"**; Title field: "Kinetic Modelling of the Maillard Reaction between Glucose and Glycine"; Author field: "Martinus A J S van Boekel"; Creator: Elsevier; distilled 2013 |
| **The year in the file name is the BOOK's, not the work's** | The abstract's own two references are Labuza & Baisier **1992** and van Boekel **1996**, and the neighbouring abstract's is Shinada et al. **1995**. This is the proceedings volume of the **6th International Symposium on the Maillard Reaction** (the conference met in 1997; the volume was first published by the Royal Society of Chemistry and reissued by Woodhead in **2005**, which is the date the PDF metadata carries). **The poster is from the late 1990s. It is not a 2005 document, and it is not the abstract of the 2005 paper — it is the announcement of the work that became it, written years earlier.** |
| No DOI | none printed; none in the metadata |
| Companions on disk | `martins2005_extraction.md` (**the paper this abstract anticipates**: Martins & van Boekel 2005, Food Chem. 90(1-2):257-269), `martins2005b_extraction.md` (the pH and initial-concentration companion), `martins2003_extraction.md` / `martins2003b_extraction.md` / `martins2003c_extraction.md` (the Amadori-degradation papers), `knol2005_extraction.md` (the same laboratory's acrylamide network, which cites this machinery) |

## 1. Why it matters

`src/kinetic_core/parameters.py`'s `MARTINS_M4` block — ten constants, `_MARTINS_SOURCE` =
"Martins & van Boekel 2005, Food Chem 90(1-2):257-269 Table 2 (= Martins 2003 thesis Table 5.2,
p. 122), model M4" — is the trunk of this engine. Its ten `_martins(...)` rows carry
`ph_of_measurement=6.8` and `temperature_range_c=(80.0, 120.0)`, and `T_REF_K = 373.15` (100 °C).
`results/validation/kinetic_core_b1_fit_report.json` is the fit built on them, and the browning
hold-out that comes out of it is the model's one out-of-sample success.

**This abstract is the same laboratory, the same author, the same system and the same programme,
written before the model existed.** Its value is therefore not numerical — it has no numbers to
give — but evidential, in three ways.

1. **It fixes the provenance.** The trunk's chemistry is not an assemblage: it is one worker's
   deliberate programme, announced as such. "Many papers pay attention only to one or two reaction
   steps or reaction products. This limits the possible kinetic analysis of the complete reaction
   network. We therefore attempted to follow the various reaction steps as closely as possible."
   The multiresponse method is cited to van Boekel 1996 (Neth. Milk Dairy J. 50:245-266), which is
   the same citation Knol 2005 and every Gökmen paper in this corpus make.
2. **It records the design decision that gives the trunk its Amadori branch.** "The Amadori
   product ... was synthesized and also heated **in the absence of sugar**, in order to be able to
   study the degradation of the Amadori product without interference of the sugar." That separate
   experiment is what identifies Martins' steps 4, 6 and 7 — the three Amadori exits the trunk
   carries as `k_ama_tdg`, `k_ama_mgo` and `k_ama_odg` — independently of the condensation that
   feeds them. It is the reason those three constants are better determined than the ones that
   depend on an unmeasured pool, and it is worth having on record next to
   `results/validation/kinetic_core_b21_prereg.md` §6, whose "unforeseen finding" is that a new
   drain on the Amadori compound worsens exactly this measured series by 24 % in half sum of
   squares. **Martins' Amadori series is not an inference from a fitted pool; it is a directly
   measured degradation of a synthesised compound.** That raises, rather than lowers, the weight the
   B21b joint fit should give it.
3. **It flags three condition mismatches against the shipped block** (section 4), one of which —
   the temperature floor — matters for how far down the trunk's barriers may be trusted.

What this abstract does **not** contain: any rate constant, any activation energy, any activation
enthalpy or entropy, any concentration, any table, any figure, any reaction scheme, any temperature
except the range, any time except "hours to days", and any result beyond a list of identified
products.

## 2. Methods as they matter to a model

Everything the abstract says about method, in its own words and in full — there is not much of it,
and all of it is quoted or closely paraphrased because the whole page is about 350 words.

- **System.** "The Maillard reaction between glucose and glycine." **No concentration is given for
  either reactant, and none is derivable.** (The 2005 paper's 0.2 mol/L each is *not* in this
  abstract.)
- **Medium and pH.** "Experiments were done in an aqueous system, **buffered to pH 7**." The buffer
  is not named and its molarity is not given. Nothing is said about whether the pH held; but the
  results paragraph says it did not — "Formation of organic acids was especially marked at the
  higher temperatures, and causes, of course, a pH decrease. As the pH strongly affects the Maillard
  reaction, as well as isomerization, this pH decrease complicates the reactions even further."
- **Temperature.** "over the temperature range **40-120 °C**."
- **Time.** "with reaction times varying from **hours to days**, depending on the temperature."
- **Vessel.** **Not stated.** No tube, no headspace, no oil bath, no stirring.
- **Quantification, response by response** (this is the abstract's most informative sentence):
  - **sugars and organic acids by HPLC**;
  - **glycine and the Amadori product, "fructosyl-lysine", by amino acid analyser**;
  - **brown colour formation by spectrometry**;
  - **fluorescent compounds by fluorometer**;
  - **"pH was measured as a function of reaction time."**
- **The separate Amadori experiment.** "The Amadori product, fructosyl-lysine, was synthesized and
  also heated **in the absence of sugar**, in order to be able to study the degradation of the
  Amadori product without interference of the sugar."
- **Products identified at the time of writing.** "The reaction products identified so far are,
  besides glucose and glycine, **fructose, fructosyl-lysine, acetic acid and formic acid**." That is
  four products. **No deoxyosone, no methylglyoxal, no melanoidin as a quantified species** — the
  3-deoxyglucosone, 1-deoxyglucosone and methylglyoxal responses that carry Martins' steps 4-8 are
  **not yet in this list**, which is the clearest sign of how early the abstract is.
- **Status of the model.** Verbatim, the whole final paragraph: "Research is underway to propose a
  kinetic model based on the main reaction products identified and the mass balance. The model
  proposed will be tested, and, if necessary, adjusted by multiresponse modelling, which also yields
  the relevant kinetic parameters (**rate constants, activation enthalpies and entropies** for every
  reaction step)."
- **Reference temperature of any fitted constant.** **There is no fitted constant, so there is no
  reference temperature.** The house rule is satisfied vacuously.
- **References, both of them, as printed.**
  1. Labuza, T. P.; Baisier, W. M. (1992). In *Physical Chemistry of Foods*, Schwartzberg, H. G.,
     Hartel, R. W., IFT Basic Symposium Series 7, Marcel Dekker, New York, pp. 595-649.
  2. Van Boekel, M. A. J. S. (1996). *Neth. Milk Dairy J.*, **50**, 245-266.

## 3. Tables re-typed

**There are none.** The abstract has no table, no figure, no scheme and no numeric result. The only
numbers on the page belonging to this abstract are the pH (7), the temperature range (40-120 °C),
the postal address, the page number (451) and the two reference citations. All of them are typed in
section 2 or section 0.

For completeness, and because it shares the page: the second abstract on p. 451 is "Improvement of
Diabetes Mellitus Complications by Dietary Antioxidants" by Nguyen van Chuyen, Jimaima Veisikiaki
Jale, Keiko Shinada, Nami Kemmotsu and Hanae Arai (Department of Food and Nutrition, Japan Women's
University, Tokyo 112, Japan). It is an in-vivo rat study on vitamin E, β-carotene, astaxanthin and
catechin against glycation endpoints, feeding for 2, 4 or 7 months; it contains no kinetics and no
concentration either, and it is unrelated to this repository's chemistry except that it measures
pentosidine and ketoamine. It is recorded here only so that a future reader who opens
`vanboekel2005.pdf` and finds two abstracts knows both were seen.

## 4. Kinetic numbers the repository can use

**None. There is not one kinetic number in this document.**

**Registry mapping (`data/keys/compounds.yml`).** Of the compounds this abstract names — glucose,
glycine, fructose, the Amadori compound, acetic acid, formic acid — **not one is in the registry**,
which is a product/marker list and carries no Maillard reactants or short-chain acids.

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| — | pH of the system | **7**, "buffered", buffer unnamed | pH unit | aqueous, 40-120 °C | — | Methods sentence, p. 451 | **level_only** (a condition, not a measurement of the system's behaviour) |
| — | temperature range | **40 to 120** | °C | aqueous, buffered pH 7 | — | p. 451 | condition |
| — | reaction times | "hours to days, depending on the temperature" | — | " | — | p. 451 | condition (no number) |
| — | products identified at the time of writing | glucose, glycine, **fructose, "fructosyl-lysine", acetic acid, formic acid** | — | " | — | p. 451 | **measured presence only — no amount** |
| — | direction: organic acid formation vs temperature | "especially marked at the higher temperatures" | — | " | — | p. 451 | within_study_ratio (a direction, no number) |
| — | direction: pH falls as acids form | stated, with the consequence that "the pH decrease complicates the reactions even further" | — | " | — | p. 451 | within_study_ratio |
| — | rate constants, barriers, enthalpies, entropies | **not yet obtained** — "Research is underway" | — | — | — | p. 451 | **absent by the author's own statement** |

### Does any number here differ from the `MARTINS_M4` block the engine ships? Three do.

The shipped block's conditions string is: *"glucose 200 mmol/L + glycine 200 mmol/L, 0.1 mol/L
phosphate, pH 6.8 initial (uncontrolled during heating), aqueous, screw-capped glass tubes, oil
bath, 80-120 C"*, with `ph_of_measurement=6.8` and `temperature_range_c=(80.0, 120.0)` on every one
of the ten rows, and `T_REF_K = 373.15` (100 °C). Set against this abstract:

| item | this abstract | `MARTINS_M4` / `martins2005_extraction.md` | verdict |
|---|---|---|---|
| **pH** | **"buffered to pH 7"** | **6.8 initial**, 0.1 mol/L phosphate, uncontrolled, falling to ≈ 5.5 in 4 h at 100 °C | **DIFFERS by 0.2 pH units.** Almost certainly the abstract rounding 6.8 to "7" for a poster; but it is a difference, it is not resolvable from either document alone, and it is recorded rather than smoothed away. **The shipped value should stay 6.8**, which is the one the full paper prints and the one the fit was made at. |
| **temperature range** | **40-120 °C** | **80, 90, 100, 110, 120 °C** (five isotherms) | **DIFFERS at the floor by 40 °C.** The abstract announces a window twice as wide at the bottom as the one that was published. Either the 40-70 °C runs were made and did not enter the M4 fit, or the plan changed. **Consequence for the engine:** the shipped `temperature_range_c=(80.0, 120.0)` is right for the fit, and the existence of unpublished 40-70 °C data is a reason to *ask*, not a licence to extrapolate. |
| **name of the Amadori compound** | **"fructosyl-lysine", printed three times** | N-(1-deoxy-D-fructos-1-yl)-**glycine** (DFG) — the transformation string in `parameters.py` is "Glc + Gly -> Schiff base (condensation)" into `AMA`, and `martins2005_extraction.md` names DFG | **THE ABSTRACT IS WRONG, and it is wrong in print, not in the OCR.** Verified on the page image at 200 dpi: the words are "fructosyl-lysine" in all three places. The Amadori compound of glucose + **glycine** is fructosyl-glycine; fructosyl-lysine is the Amadori compound of glucose + **lysine** (the milk / protein compound, the one `hamzalioglu2026_extraction.md` deals with as lactulosyl-lysine's sibling). This is a slip in a poster abstract about a glucose/glycine system whose title says glycine. Flags 2. |
| **reactant concentrations** | **not given** | 200 mmol/L each | the abstract simply does not say; no conflict |
| **buffer** | "buffered", not named | 0.1 mol/L phosphate | the abstract does not say; no conflict |
| **vessel** | not stated | 10 mL screw-capped Schott tubes, oil bath | the abstract does not say; no conflict |
| **responses** | sugars, organic acids, glycine, Amadori, **brown colour**, **fluorescent compounds**, pH | the M4 fit's responses are glucose, fructose, the Amadori compound, formic acid, acetic acid, methylglyoxal, the two deoxyosones and melanoidin | **the abstract lists a FLUORESCENCE response the published model does not use**, and lists neither deoxyosone nor methylglyoxal. Flags 3. |
| **temperature-dependence formalism** | "activation **enthalpies and entropies**" (an Eyring reading) | activation **energies** in a reparameterised Arrhenius, k = X·exp(−Y·Ea) with X at T_av | the plan and the publication differ in formalism; the shipped `ea_kj_mol` values are Arrhenius activation energies and must not be described as enthalpies. Flags 4. |
| **any rate constant** | **none** | ten, with 95 % HPD on both k and Ea | nothing to compare |

**So: no number in this abstract contradicts a shipped constant, because it contains no constants.
Two of its three condition statements differ from the shipped conditions (pH 7 vs 6.8; 40-120 vs
80-120 °C), and its third — the identity of the Amadori compound — is an outright error in the
source.** Nothing in `parameters.py` should change on the strength of this document. It should be
cited, if at all, only as provenance for the programme, never as a source of a value.

## 5. Flags

1. **The file name's year is the reprint's, not the work's.** `vanboekel2005.pdf` is page 451 of
   *The Maillard Reaction in Foods and Medicine*, whose 2005 imprint the PDF metadata records; the
   abstract's own references stop at 1996, and its companion abstract's at 1995. **Anything that
   cites this as "van Boekel 2005" implies a 2005 document and will mislead.** The safe citation is
   "van Boekel, poster abstract, *The Maillard Reaction in Foods and Medicine*, p. 451". It is **not**
   a short version of Martins & van Boekel 2005 (Food Chem. 90:257-269), which is a different
   document with a different first author, and which is what `_MARTINS_SOURCE` correctly names.
2. **"Fructosyl-lysine" for a glucose/glycine Amadori compound is an error printed three times, and
   it was verified on the page image, so it is not an OCR artefact.** Do not propagate it. The
   compound is N-(1-deoxy-D-fructos-1-yl)-glycine. The risk is specific and live: the repository
   also carries a genuine fructosyl-lysine lane (wave B20's `FLP`, from Nguyen 2016 and Berk 2021,
   and Hamzalıoğlu 2026's lactulosyl-lysine, which is the constant B21 transferred to water). **A
   text search for "fructosyl-lysine" across the corpus will match this abstract and mean something
   entirely different.** That is worth a note wherever the two lanes are indexed together.
3. **The abstract lists a response the published model never used: fluorescent compounds by
   fluorometer.** Fluorescence is an early-browning readout independent of A420/A470, and no
   fluorescence datum from this laboratory is anywhere in this corpus. Since the trunk's browning
   hold-out is its one out-of-sample success and rests on a single absorbance response, **a second,
   independent browning readout on the same pots would be a genuinely new test.** It is the single
   most valuable thing this one page tells us to ask for.
4. **The plan was Eyring, the publication was Arrhenius.** "Activation enthalpies and entropies"
   here; activation energies in Table 2 there. The shipped `ea_kj_mol` fields are activation
   energies. If an activation entropy for any of the ten steps exists anywhere, it is in the 2003
   thesis and not in the 2005 paper; it is not on disk.
5. **The product list is four compounds long and does not yet include the deoxyosones or
   methylglyoxal.** Fructose, the Amadori compound, acetic acid and formic acid only. Martins'
   published steps 4 through 8 — the ones the trunk depends on for 3-deoxyglucosone,
   1-deoxyglucosone and methylglyoxal — rest on responses that did not exist when this was written.
   That is a dating aid, not a criticism; it places the poster before the 2003 thesis work.
6. **What this document does not contain**: any rate constant; any barrier; any enthalpy or entropy;
   any concentration of anything; any reactant loading; any buffer identity or molarity; any vessel
   description; any table; any figure; any reaction scheme; any temperature other than the range;
   any time other than "hours to days"; any statement of replication; and any result that could be
   scored.
7. **What to request**: (i) whether the **40 to 70 °C** runs the abstract announces were made, and
   if so where their data are — that window is entirely absent from the corpus and would be the
   only check on how far Martins' barriers may be extrapolated downward, which is the assumption
   every low-temperature answer from this engine rests on; (ii) the **fluorescence** data (Flags 3);
   (iii) confirmation that the buffer was the 0.1 mol/L phosphate at pH 6.8 the 2005 paper reports,
   so that the "pH 7" here can be closed out as rounding.
8. **Registry gaps against `data/keys/compounds.yml`**: none of glucose, glycine, fructose, the
   Amadori compound, acetic acid or formic acid is keyed. This is the same gap
   `knol2005_extraction.md` §5 item 12 records; the registry carries products and markers, not
   Maillard reactants or C1-C2 acids, and a benchmark built on the trunk's own source paper would
   need at least glucose, glycine and the two acids keyed.
