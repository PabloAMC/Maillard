# Tutorial: three worked examples, and how to read the output

You do not need to read the source to use this model. You do need to read this page, because
the model is unusually opinionated about what it will and will not tell you, and most of the
value is in the parts where it refuses. Install and the command reference are in the
[quick start](guides/QUICKSTART.md); the chemistry and how well the model does are in the
[introduction](guides/INTRODUCTION.md).

**One line, before anything else:** compare two formulations and read the **ratio**. Never quote
an absolute ppb number as a specification. Every ratio is printed with a reliability grade for the
axis the two arms differ on; that grade, not the ratio's size, says whether to act on it.

---

## 1. What a spec looks like

A spec is YAML (or JSON; the same loader reads both). Two arms, `a:` and `b:`, for a comparison;
one flat block for a single prediction. This is `docs/examples/compare_ribose_vs_glucose.yml`:

```yaml
a:
  name: cysteine_ribose
  precursors:
    L-Cysteine: 10.0      # MILLIMOLAR in the reacting phase. Not grams, not %w/w.
    D-Ribose: 10.0
  temp_C: 140.0           # measured product temperature, not oven set point
  time_min: 30.0          # length of the hold at temp_C
  ph: 5.0
  aw: 0.98
  protein_type: free      # 'free' / 'water' = aqueous model system
b:
  name: cysteine_glucose
  precursors:
    L-Cysteine: 10.0
    D-Glucose: 10.0
  temp_C: 140.0
  time_min: 30.0
  ph: 5.0
  aw: 0.98
  protein_type: free
```

Optional keys worth knowing:

| key | what it does |
| --- | --- |
| `targets:` | ask for named compounds instead of the lane's defaults. Asking for something the model cannot represent produces a **named refusal with its reason**, which is often the most useful output you can get. |
| `matrix:` | a matrix on file. Two kinds: an odour-threshold matrix (`gelatin_3pct`, `skim_milk`) selects that matrix's **measured** thresholds where they exist and says "no measured threshold" where they do not; a protein matrix (`blg`, `soy_isolate`, `pea_isolate`) together with `protein_g_per_l` charges the protein's reactive sites, so the thiols meet its disulfide pool and the aldehydes its amine pool. Nothing is ever borrowed from another matrix. |
| `protein_sites:` | your own isolate's site densities in mmol per gram (free thiol, disulfide, amine), for a protein not on file. The quick start shows both forms. |
| `measured_matrix_ratios:` | your own measured matrix-to-water ratios, per compound. Turns on the report's **residual decomposition**: how much of *your* measured shift the model's named terms explain, and how much is unexplained. |

**Units are not negotiable.** `precursors` is mM. A protein isolate is not a precursor and has no
molar basis: state it as a loading in grams per litre with a matrix on file or its own sites, and
the model charges what it can (the sites) and says out loud what it cannot (the isolate's own
sugars, free amino acids and lipid volatiles, which it does not know).

---

## 2. Three worked examples

### Example 1: the decision. Does ribose beat glucose?

```bash
python scripts/maillard.py compare docs/examples/compare_ribose_vs_glucose.yml \
    --report /tmp/ribose_vs_glucose.html
```

Output, abridged (the numbers are the model's on 9 September 2026; run it to see today's):

```text
  COMPARE [KINETIC CORE]   A = cysteine_ribose   vs   B = cysteine_glucose

  arm A:
  ENVELOPE: IN_ENVELOPE_EXTRAPOLATED   lane: sulfur
  mapped precursors: Cys=10 mM, PENT=10 mM
    ~ declared extrapolation -- no buffer was declared for this system, so the pH TRAJECTORY is
    ~ EXTRAPOLATED: it is computed from water autoprotolysis and the charged solutes alone. ...

  arm B:
  ENVELOPE: IN_ENVELOPE_EXTRAPOLATED   lane: sulfur
  mapped precursors: Cys=10 mM, Glc=10 mM
    ~ declared extrapolation -- HEXOSE ENTRY UNIDENTIFIED (FFT, MFT): the only route from a hexose
    ~ to these thiols is the C2+C3 fragmentation entry, whose rate constants no primary
    ~ measurement identifies ... the ordering 'pentose above hexose' is the structural claim the
    ~ model does support.

  Axes this comparison moves: sugar_identity
  Governing reliability:      do-not-use (4/10)   [directional panel, independent claims]

  1 of 6 ratios resolve above the same-sample dispersion band (4.8x)
  ...of which 4 are UNDEFINED (one arm at exactly zero, or on a route the engine declares
  unidentified) and resolve nothing: see the 'resolved' column, not this count.

  compound                                        A/B direction      resolved
  2-acetylthiazole                            0.0251x higher_in_cysteine_glucose yes
  2-furfurylthiol (FFT)                  B unidentified undefined      n/a -- undefined
  2-methyl-3-furanthiol (MFT)            B unidentified undefined      n/a -- undefined
  bis(2-methyl-3-furyl) disulfide        B unidentified undefined      n/a -- undefined
  furfural                                     0.923x higher_in_cysteine_glucose no -- inside band
  methanethiol                                 A only undefined      n/a -- undefined
```

**How to read it.** Only the `resolved` rows are claims, and here there is one. The two thiols
and their disulfide are **undefined**, not "higher in ribose": the glucose arm's only route to
them runs through a step no measurement identifies, so its number is a floor artefact and no
ratio is claimed. What the model does support is the ordering, pentose above hexose, and it says
so in the arm's declaration rather than printing a ratio of a million. Furfural's 0.92× sits
inside the same-sample dispersion band and is **NOT RESOLVED**: not "a small effect", but "a
difference this analytical method cannot see". The report styles those rows in grey with the
ratio struck through, so you cannot quote one by accident. `A only` means arm B is at exactly
zero, so there is no ratio at all.

The reliability line is the other half. The axis this comparison moves is sugar identity, and the
model's direction sense on that axis is graded from the independent claims of the directional
panel; today the grade is *do-not-use*, and the line names the claims it missed. A ratio you can
act on is one with `resolved: yes` **and** a grade of *caution* or better. Neither arm declared a
buffer, so both pH trajectories are extrapolated; because both arms share that extrapolation, the
ratio is more trustworthy than either absolute, which is the whole argument for reading ratios.

### Example 2: the profile. What does a cysteine and ribose reaction flavour smell of?

```bash
python scripts/maillard.py predict docs/examples/compare_ribose_vs_glucose.yml \
    --system a --report /tmp/cys_ribose.html
```

```text
  compound                                  ug/L (= ppb in water)             OAV
  furfural                                                 135.77          0.0453
  2-methyl-3-furanthiol (MFT)                               76.27        1.53e+04
  2-furfurylthiol (FFT)                                     57.67        9.61e+03
  2-acetylthiazole                                           2.10            0.21
  bis(2-methyl-3-furyl) disulfide                            0.12             376
  methanethiol                                                  0    no threshold
```

**How to read it.** Read the OAV column, not the µg/L column. Furfural is the most abundant
compound and contributes essentially nothing (OAV 0.05); the two thiols are less abundant and sit
four orders of magnitude above their thresholds. Abundance is not aroma. In the HTML report each
number carries its interval, roughly a fiftyfold band, which is the honest width of an absolute
here, and the OAV chart plots them on a log axis with whiskers. The disulfide is charted at its
potency-weighted value, because it is far more potent than its own monomer: mass lost to
dimerisation is not aroma lost. Methanethiol has no measured threshold in water anywhere in the
corpus, so it is reported as having none rather than being given a borrowed one.

Before you believe any of these absolutes, read the model card in the README: out of sample, a
small minority of panel rows land within threefold. The card is generated from the artifacts, so
this page does not repeat its numbers.

### Example 3: the refusal, which is also an answer

```yaml
# refusal_demo.yml
name: cys_ribose_asking_for_hmf
precursors: {L-Cysteine: 10.0, D-Ribose: 10.0}
targets: ["2-methyl-3-furanthiol (MFT)", "HMF", "2-pentylfuran"]
temp_C: 140.0
time_min: 30.0
ph: 5.0
aw: 0.98
```

```bash
python scripts/maillard.py predict refusal_demo.yml --report /tmp/refused.html
```

```text
  ENVELOPE: OUT_OF_ENVELOPE   lane: sulfur
  mapped precursors: Cys=10 mM, PENT=10 mM
    ! REFUSED -- UNREPRESENTED TARGETS: 2-pentylfuran -- The lipid lane exists, but 2-pentylfuran
    ! is NOT in Frankel 1989's six-product slate and no branch fraction for the linoleate ->
    ! alkylfuran route is measured anywhere in the fit corpus. The retired screening lane's
    ! shipped 0.08 had no source. Refused rather than invented.
    ~ declared extrapolation -- no buffer was declared for this system, so the pH TRAJECTORY is
    ~ EXTRAPOLATED ...
    ~ declared extrapolation -- 5-HMF: the two formation limbs are ingested WHOLE from Kocadagli &
    ~ Gokmen 2016's AMINE-FREE amorphous glucose melt at 160-200 C. This program runs at 140 C in
    ~ an aqueous or matrix system, so both the temperature and the physical state are
    ~ extrapolations. ...
    ~ declared extrapolation -- 5-HMF: THE MODEL HAS NO VALIDATED SINK AT COOKING TEMPERATURE. ...
    ~ The furanic extraction dossier's declared gap G2: the 50-150 C window is empty. EXPECT HMF
    ~ TO BE OVER-PREDICTED.
    ~ declared extrapolation -- 5-HMF + cysteine: the sink constant is HELD at its 50 C value for
    ~ this whole program. Holding it UNDER-states the sink; extrapolating it is a named prohibited
    ~ derivation, and the direction is stated rather than chosen for convenience.

  NO NUMBER IS EMITTED. The core declined this request above.
```

and on stderr the same refusal, closing with:

```text
  A refusal is an output, not a failure. Run `python scripts/maillard.py explain <compound>` to
  see what the core does carry for a compound, and why.
```

**This is the feature, not the failure.** The alternative is a plausible-looking number with
nothing behind it, and every documented accuracy defect in this repository began as a number that
should not have existed. The refusal tells you what would have to be measured for the answer to
exist, which is a research plan, not an error message.

**Watch what happened to HMF, because it is the other half of the lesson.** Until the furanic
channels were added, this same spec refused *two* targets, and the HMF refusal read "the
hexose-dehydration route that forms it was never parameterised". That step parameterised it, so
HMF is now an answerable species and only 2-pentylfuran refuses. **The refusal did not become a
silent pass; it became four declared extrapolations**, one of which states the expected direction
of the error out loud. A refusal is what the model says when it has no route; a declared
extrapolation is what it says when it has a route it does not trust at your conditions. Neither
is a number you should ship. The HTML report renders each refusal as its own card, beside the
list of what the model *can* be asked and every compound it deliberately refuses, with the reason.

---

## 3. Reading the outputs

### Intervals

Every absolute the core emits is an **interval**, never a point. The floor width is set by two
measured facts:

- **Same-sample dispersion of headspace GC-MS, 10 to 23 times.** Two papers measuring the *same
  samples* disagree by that much. It is a calibration fact, not a fitted error.
- **Half a decade on the air-to-water partition constant.** The literature spread on hexanal's
  partition constant alone is nearly tenfold; the ruling is to *keep* the band, not to pick a
  value.

Added in quadrature, that is roughly a fiftyfold band before any model error at all. A lipid-lane
compound has a wider one, because its rate is an assumption; the report's declared-assumptions
section lists exactly which assumptions widened your run and by how much.

**A compound with no interval is weaker evidence than one with a wide interval, not stronger.**

### NOT RESOLVED

A ratio inside the same-sample dispersion band (about 0.2× to 4.8×) is reported NOT RESOLVED. The
model is not saying the two arms are the same; it is saying this method cannot tell them apart.
Do not report it as a null result.

### Refusals and undefined ratios

Five kinds, none of which emits a number:

| refusal | what it means |
| --- | --- |
| **unmapped precursor** | you named something the core cannot charge: an intact protein as a precursor, a flour, an isolate in mM |
| **unrepresented target** | you asked for a compound the core cannot name. Today the list is 1-hexanol, 2-pentylfuran, propanal, 2-nonenal, HEMF (homofuraneol) and 2,5-dimethyl-4-hydroxy-3(2H)-thiophenone; the engine's `UNREPRESENTED_COMPOUNDS` carries each one's reason. HMF, DMHF and furaneol left this list when the furanic channels were added and now answer. |
| **lane conflict** | your request needs two Maillard lanes at once; they do not compose, because that would spend the same cysteine twice |
| **missing sulfur source or asparagine** | the lane was selected but the charge cannot supply the atom the product is made of |
| **unidentified route** (a comparison) | one arm reaches the compound only through a step no measurement identifies, as the glucose arm reaches the thiols in Example 1; the ratio is printed as *undefined* with the arm named, never as a magnitude |

Run `python scripts/maillard.py explain <compound>` to see what the model does have, and why.

### Declared extrapolations

Distinct from a refusal: the model **answered**, and is telling you the answer sits outside what
its parameters license. The commonest are "no buffer declared, so the pH trajectory is
extrapolated" and, on the lipid lane, "the rate anchor was measured at 25 °C and this program
peaks at 150 °C". They are not footnotes; they are the reason the interval is that wide.

---

## 4. Where compounds come from: `explain`

```bash
python scripts/maillard.py explain MFT
python scripts/maillard.py explain HMF             # two routes, one of them pinned
python scripts/maillard.py explain 2-pentylfuran   # a refusal, with its reason and the route the literature draws
```

It prints every route the model has to that compound, the **evidence class** of each step, and
the literature anchors those steps rest on:

| class | meaning |
| --- | --- |
| `measured` | a rate constant or activation energy printed in a paper |
| `fitted` | estimated here by least squares on declared fit rows |
| `derived` | computed from another constant by a stated relation |
| `pinned` | held at a value nothing measured, including held at zero |

`explain HMF` is worth running for the contrast: two formation routes, one measured and one
pinned, and the pinned one keeps the authors' own zero activation energy, quoted verbatim from
their footnote, because no defensible value for that edge exists in any paper of the cluster.
`explain 2-pentylfuran` prints no route the model *has*, the reason it is refused, the list of
compounds the model *can* explain, and the cited reaction rule by which the hypothesis layer
reaches it from the linoleate hydroperoxides: the refusal reads "no rate, not no route".

Nothing on that page is new data; every line is read from a frozen registry at run time, so it
moves the day the model moves. If you are deciding whether to believe a prediction, this is the
command that answers it.

---

## 5. The network map

```bash
python scripts/generators/generate_network_map.py --all-examples
# -> docs/assets/network_map.html
```

All four lanes, drawn from the live code: nodes are species, edges are reactions coloured by
evidence class, and every edge's tooltip shows its source anchor and validity window. Dim
everything except `pinned` in the legend to see how much of the network rests on constants nothing
measured.

**Flux mode** additionally sets edge width by the time-integrated flux through each step for a
given process, so you can see where *your* chemistry actually went. Three example processes ship
with it: a 145 °C cysteine and ribose reaction flavour, a three-zone pea-isolate extrusion (sulfur
and lipid co-integrated), and a 180 °C acrylamide fry. The same network looks completely different
in each. Point it at your own process with `--spec my_process.json --out my_map.html`; flux specs
accept `thermal_segments: [[duration_min, temp_C], ...]`, so a multi-zone process can be written as
the programme it physically is.

---

## 6. Your own data

`score` scores your pots the way the panel is scored and refits nothing; `calibrate` writes a
per-laboratory file from the same document (your levels set a response factor per compound, your
contrasts move the few rate constants they can identify, the hold-out is scored before and after),
which `--calibration` applies to the other verbs while the shipped model stays untouched. Both are
walked through, with a real laboratory's four-temperature ladder as the example, in the quick
start's [score](guides/QUICKSTART.md#your-own-measurements-maillard-score) and
[calibrate](guides/QUICKSTART.md#calibrate-to-your-laboratory-maillard-calibrate) sections. What a
calibration cannot do: fix a structural miss. If the model gets the sign of a temperature trend
wrong, a response factor moves the level and the card will show the trend still wrong.

---

## 7. What this model is good for, and what it is not

### Use it for

- **Choosing between two formulations.** Ratios cancel the shared systematic error. The
  directional panel's current score on strictly independent claims is the first line of
  `python scripts/maillard.py --help` and of the README model card; this page does not repeat it.
- **Sugar swaps and cysteine present against absent**, read with the per-axis grade the `compare`
  verb prints beside every ratio. That column comes from the same panel; read it rather than a
  number quoted here.
- **Deciding what to measure next.** `rank` orders candidate measurements by how much each would
  move the model, and `wishlist` says what each would unlock. Every ranked row is a place the
  model is *measurably wrong*, so this is the one claim type that does not depend on the model
  being right.
- **Auditing a mechanism.** `explain` and the network map show what the model believes and on what
  evidence, including where it believes something because nothing measured it.

### Do not use it for

- **Any absolute ppb as a specification.** The model card's out-of-sample line is the number to
  quote, and it is small.
- **A direction on an axis whose grade is *do-not-use*.** The grade is printed with every ratio.
  A lane with no term for the axis you moved (the lipid lane has neither a pH nor a
  water-activity term) refuses the comparison rather than answering with identical arms.
- **Sulfur absolutes at long times or low temperatures.** The model removes the thiols too fast
  at 100 °C and at 140 °C; the introduction's section 7 is about exactly this, and both structures
  tried so far for the missing sink were refused. The directions at 145 °C are a separate, better
  question.
- **Anything that happens after the cook.** The model has no storage clock. The one post-cook
  process it knows about is the covalent aldehyde-to-protein sink, and it knows about it well
  enough to tell you it is the wrong tool: that channel's activation energy was measured at 15 to
  23 kJ/mol, so it removes a fraction of a percent of an aldehyde during any real thermal step,
  but over weeks at ambient in a high-protein matrix it is real and sizeable. If your question is
  about shelf life rather than the cook, this model answers the wrong question, and it will not
  warn you, because every number it prints is an end-of-cook number.
- **Anything the model refuses.** Do not work around a refusal by substituting a nearby compound.
  Norfuraneol is not DMHF; that substitution is exactly what the refusal exists to prevent.

### The single most useful habit

Run the comparison, read the resolved ratios and their grade, ignore the absolutes, and book the
experiment that `rank` puts at the top. The model is a well-instrumented way of being wrong in a
located, quantified way, and that is worth more than a confident number.

---

*Authoritative accuracy figures live in the generated model card in the
[README](../README.md#how-well-calibrated-is-it), which is regenerated from the artifacts. If this
page and the model card ever disagree, the model card is right and this page has drifted.*
