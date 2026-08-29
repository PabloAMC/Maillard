# Using the tool — a quickstart for a flavour scientist

You do not need to read the source to use this model. You do need to read this page,
because the model is unusually opinionated about what it will and will not tell you, and
most of the value is in the parts where it refuses.

**One line, before anything else:** compare two formulations and read the **ratio**. Never
quote an absolute ppb number as a specification. Treat any pH or moisture direction as
unsupported.

---

## 0. The 60-second version

```bash
# what the tool is and how to ask it things
python scripts/maillard.py --help

# a starter spec you can edit
python scripts/maillard.py compare --template > my_experiment.yml

# the thing you actually want: A vs B, with an HTML report
python scripts/maillard.py compare my_experiment.yml --report report.html
```

Open `report.html` in any browser. It is a single file with no network dependency — mail it,
print it, open it on a plane.

---

## 1. What a spec looks like

A spec is YAML (or JSON — the same loader reads both). Two arms, `a:` and `b:`, for a
comparison; one flat block for a single prediction.

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
| `matrix:` | a food matrix (`gelatin_3pct`, `skim_milk`, …). Selects that matrix's **measured** odour thresholds where they exist and an explicit "no measured threshold" where they do not. Nothing is ever borrowed from another matrix. |
| `measured_matrix_ratios:` | your own measured matrix/water ratios, per compound. Turns on the report's **residual decomposition**: how much of *your* measured shift the model's named terms explain, and how much is unexplained. |

**Units are not negotiable.** `precursors` is mM. A protein isolate is not a precursor — it is
a *lipid carrier*, its declared charge is ignored (mM of an isolate has no molar basis), and
the model says so out loud on every run that uses one.

---

## 2. Three worked examples

### Example 1 — the decision: does ribose beat glucose?

```bash
python scripts/maillard.py compare data/cli_examples/compare_ribose_vs_glucose.yml \
    --report /tmp/ribose_vs_glucose.html
```

Expected output (abridged):

```text
================================================================================
  COMPARE [KINETIC CORE]   A = cysteine_ribose   vs   B = cysteine_glucose
================================================================================

  arm A:
  ENVELOPE: IN_ENVELOPE_EXTRAPOLATED   lane: sulfur
  mapped precursors: Cys=10 mM, PENT=10 mM
    ~ declared extrapolation -- no buffer was declared for this system, so the
    ~ pH TRAJECTORY is EXTRAPOLATED: it is computed from water autoprotolysis and
    ~ the charged solutes alone. ...

  arm B:
  ENVELOPE: IN_ENVELOPE_EXTRAPOLATED   lane: sulfur
  mapped precursors: Cys=10 mM, Glc=10 mM
    ~ declared extrapolation -- ...

  4 of 6 ratios resolve above the same-sample dispersion band (4.8x)
  ...of which 2 are UNDEFINED (one arm at exactly zero) and resolve nothing:
  see the 'resolved' column, not this count.

  compound                                        A/B direction      resolved
  ------------------------------------------------------------------------------
  2-acetylthiazole                          2.04e-05x higher_in_cysteine_glucose yes
  2-furfurylthiol (FFT)                         17.2x higher_in_cysteine_ribose  yes
  2-methyl-3-furanthiol (MFT)                    244x higher_in_cysteine_ribose  yes
  bis(2-methyl-3-furyl) disulfide              A only undefined       n/a -- undefined
  furfural                                      53.2x higher_in_cysteine_ribose  yes
  methanethiol                                 A only undefined       n/a -- undefined
```

**How to read it.** Only the `resolved` rows are claims. A ratio between 0.21× and 4.8× falls
inside the measured same-sample dispersion band and is reported **NOT RESOLVED** — that is not
"a small effect", it is "a difference this analytical method cannot see". The report styles
those rows in grey with the ratio struck through, so you cannot accidentally quote one.

`A only` means arm B is at exactly zero, so there is no ratio at all. That is a third state,
`n/a — undefined`, and it is neither a resolved claim nor an unresolved one.

The extrapolation block above the table is not boilerplate: neither arm declared a buffer, so
the sulfur lane's pH trajectory is being extrapolated in both. Because *both* arms carry the
same extrapolation, the ratio is more trustworthy than either absolute — which is the whole
argument for reading ratios.

This is the comparison the model is best at. The only axis that moves is sugar identity, which
scores 9/11 on the directional panel.

### Example 2 — the profile: what does a cysteine/ribose reaction flavour smell of?

```bash
python scripts/maillard.py predict data/cli_examples/compare_ribose_vs_glucose.yml \
    --system a --report /tmp/cys_ribose.html
```

```text
  compound                                  ug/L (= ppb in water)             OAV
  ----------------------------------------------------------------------------
  furfural                                              194.90           0.065
  2-furfurylthiol (FFT)                                 119.40        1.99e+04
  2-methyl-3-furanthiol (MFT)                            62.44        1.25e+04
  2-acetylthiazole                                        2.10            0.21
  bis(2-methyl-3-furyl) disulfide                            0               0
  methanethiol                                               0    no threshold
```

**How to read it.** Read the OAV column, not the µg/L column. Furfural is the most abundant
compound and contributes essentially nothing (OAV 0.065); FFT is 1.6× *less* abundant than
furfural and sits four orders of magnitude above its threshold. Abundance is not aroma.

In the HTML report each of those numbers carries its interval — roughly a 49× band, which is
the honest width of an absolute here — and the OAV chart plots them on a log axis with
whiskers. **The MFT dimer is charted at its potency-weighted value**, because it is ~15.6×
more potent than its own monomer: mass lost to dimerisation is *not* aroma lost.

`methanethiol` has no measured threshold in water anywhere in the corpus, so it is reported as
having none rather than being given a borrowed one.

### Example 3 — the refusal, which is also an answer

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
    ! REFUSED -- UNREPRESENTED TARGETS: 2-pentylfuran -- The B6 lipid lane exists, but
    ! 2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the
    ! linoleate -> alkylfuran route is measured anywhere in the fit corpus. The FAST lane's
    ! shipped 0.08 has no source. Refused rather than invented.
    ~ declared extrapolation -- no buffer was declared for this system, so the pH TRAJECTORY is
    ~ EXTRAPOLATED: it is computed from water autoprotolysis and the charged solutes alone. ...
    ~ declared extrapolation -- 5-HMF: the two formation limbs are ingested WHOLE from Kocadagli &
    ~ Gokmen 2016's AMINE-FREE amorphous glucose melt at 160-200 C. This program runs at 140 C in
    ~ an aqueous or matrix system, so both the temperature and the physical state are
    ~ extrapolations. ...
    ~ declared extrapolation -- 5-HMF: THE MODEL HAS NO VALIDATED SINK AT COOKING TEMPERATURE. ...
    ~ gap G2: the 50-150 C window is empty. EXPECT HMF TO BE OVER-PREDICTED.
    ~ declared extrapolation -- 5-HMF + cysteine: the sink constant is HELD at its 50 C value for
    ~ this whole program. Holding it UNDER-states the sink; extrapolating it is a named prohibited
    ~ derivation (K5a sec. 7.3), and the direction is stated rather than chosen for convenience.

  NO NUMBER IS EMITTED. The core declined this request above.
```

on stderr:

```text
OUT OF ENVELOPE -- no number is emitted. Lane resolved: sulfur.
  missing species (targets the core cannot name): 2-pentylfuran
  declared reason -- UNREPRESENTED TARGETS: 2-pentylfuran -- The B6 lipid lane exists, but
  2-pentylfuran is NOT in Frankel 1989's six-product slate and no branch fraction for the
  linoleate -> alkylfuran route is measured anywhere in the fit corpus. The FAST lane's
  shipped 0.08 has no source. Refused rather than invented.
  A refusal is an output, not a failure. Run `python scripts/maillard.py explain <compound>` ...
```

**This is the feature, not the failure.** The alternative is a plausible-looking float with
nothing behind it, and every documented accuracy defect in this repository began as a number
that should not have existed. The refusal tells you precisely what would have to be measured
for the answer to exist — which is a research plan, not an error message.

**Watch what happened to HMF here, because it is the other half of the lesson.** Until Wave B7
this same spec refused *two* targets, and the HMF refusal read "5-HMF is not a species in any core
lane. The hexose-dehydration route that forms it was never parameterised." B7 parameterised it, so
HMF is now an answerable trunk species and only 2-pentylfuran refuses. **The refusal did not become
a silent pass — it became four declared extrapolations**, one of which states the expected
direction of the error out loud ("EXPECT HMF TO BE OVER-PREDICTED"). A refusal is what the model
says when it has no route; a declared extrapolation is what it says when it has a route it does not
trust at your conditions. Neither is a number you should ship.

The HTML report renders each refusal as its own card, alongside the full list of what the
model *can* be asked and every compound it deliberately refuses with the reason.

---

## 3. Reading the outputs

### Intervals

Every absolute the core emits is an **interval**, never a point. The floor width is set by two
measured facts:

- **HS-SPME same-sample dispersion, 10–23×** — two papers measuring the *same samples* disagree
  by that much. It is a calibration fact, not a fitted error.
- **±0.5 decades on the air/water partition constant** — the literature spread on hexanal's
  K_aw alone is 9.5×; the ruling is to *carry* the band, not to pick a constant.

Added in quadrature, that is a ~49× band before any model error at all. A lipid-lane compound
carries more, because its rate is an assumption; the report's declared-assumptions section
lists exactly which assumptions widened your run and by how much.

**A compound with no interval is weaker evidence than one with a wide interval, not stronger.**

### NOT RESOLVED

A ratio inside the same-sample dispersion band (~0.21×–4.8×) is reported NOT RESOLVED. The
model is not saying the two arms are the same; it is saying this method cannot tell them
apart. Do not report it as a null result.

### Refusals

Four kinds, all of which emit no number:

| refusal | what it means |
| --- | --- |
| **unmapped precursor** | you named something the core cannot charge — an intact protein, a flour, an isolate used as a precursor |
| **unrepresented target** | you asked for a compound the core cannot name. The current list is **1-hexanol, 2-pentylfuran, propanal, 2-nonenal, HEMF / homofuraneol** and **2,5-dimethyl-4-hydroxy-3(2H)-thiophenone** — read them off `engine.UNREPRESENTED_COMPOUNDS`, which carries each one's reason. *HMF, DMHF and furaneol left this list in Wave B7 and now answer.* |
| **lane conflict** | your request needs two Maillard lanes at once; they do not compose, because that would spend the same cysteine twice |
| **missing sulfur source / asparagine** | the lane was selected but the charge cannot supply the atom the product is made of |

Run `python scripts/maillard.py explain <compound>` to see what the model does carry, and why.

### Declared extrapolations

Distinct from a refusal: the model **answered**, and is telling you the answer sits outside
what its parameters license. The commonest ones are "no buffer declared, so the pH trajectory
is extrapolated" and, on the lipid lane, "the rate anchor was measured at 25 °C and this
program peaks at 150 °C". They are not footnotes; they are the reason the interval is that
wide.

---

## 4. Where compounds come from: `explain`

```bash
python scripts/maillard.py explain MFT
python scripts/maillard.py explain hexanal
python scripts/maillard.py explain HMF             # answers since B7 -- 2 routes, 1 of them `pinned`
python scripts/maillard.py explain 2-pentylfuran   # a refusal, with its declared reason
```

`explain HMF` is worth running for the contrast: it resolves to `lane: trunk`, two formation
routes, **measured 1 / pinned 1** — and the pinned one carries the authors' own `Ea = 0`, quoted
verbatim from their footnote, because no defensible activation energy for that edge exists in any
paper of the cluster. `explain 2-pentylfuran` prints no routes at all, the reason it is refused,
and the full list of compounds the model *can* explain.

Prints every route the model has to that compound, the **evidence class** of each step, and
the literature anchors those steps rest on:

| class | meaning |
| --- | --- |
| `measured` | a rate constant or activation energy printed in a paper |
| `fitted` | estimated here by least squares on declared FIT rows |
| `derived` | computed from another constant by a stated relation |
| `pinned` | held at a value nothing measured — including held at zero |

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
evidence class, and every edge's tooltip carries its source anchor and validity window. Dim
everything except `pinned` in the legend to see how much of the network rests on constants
nothing measured.

**Flux mode** additionally sets edge width by the time-integrated absolute flux through each
step for a given process, so you can see where *your* chemistry actually went. Three example
processes ship with it: a 145 °C cysteine/ribose reaction flavour, a three-zone pea-isolate
extrusion (sulfur and lipid co-integrated), and a 180 °C acrylamide fry. The same network looks
completely different in each.

Point it at your own process:

```bash
python scripts/generators/generate_network_map.py --spec my_process.json --out my_map.html
```

Flux specs accept `thermal_segments: [[duration_min, temp_C], ...]` so a multi-zone process can
be written as the program it physically is.

---

## 6. What this model is good for — and what it is not

### Use it for

- **Choosing between two formulations.** Ratios cancel the shared systematic error; the
  directional panel scores 24/36 on strictly independent claims, and 18/23 once pH and water
  activity are set aside.
- **Sugar swaps** (9/11 on the panel) and **cysteine present vs absent** (4/4). These are the
  model's two strongest axes.
- **Deciding what to measure next.** `rank-experiments` ranks candidate measurements by how
  much each would move the model. Every ranked row is a place the model is *measurably wrong*,
  so this is the one claim type that does not depend on the model being right.
- **Auditing a mechanism.** `explain` and the network map show what the model believes and on
  what evidence, including where it believes something because nothing measured it.

### Do not use it for

- **Any absolute ppb as a specification.** Out of sample: 6.04× median fold error on the
  free-precursor hold-out (worst 52.6×), 67–94× on the matrix lane, and 1 of 5 genuine
  extrapolation rows inside the 90 % CI.
- **pH direction** (6/10 on the panel — at or near chance) or **water activity direction**
  (0/3). The trunk and acrylamide lanes carry no pH term and no lane carries an a_w term at
  all; your value is recorded and ignored.
- **Sulfur absolutes.** The sulfur branch has 8 primary-source-verified literature anchors and
  the model fails every one of them. The *directions* are a separate question, and the
  temperature direction is also wrong.
- **The acrylamide lane's time shape.** It is inverted against measurement — a source measures
  acrylamide rising 28 → 1459 ppb between 10 and 30 min where the core predicts it falling.
- **Anything that happens after the cook.** The model has no storage clock. The one
  post-cook process it knows about is the covalent aldehyde–protein sink, and it knows about
  it well enough to tell you it is the wrong tool: that channel's activation energy was
  measured in 2026 at **15–23 kJ/mol**, so it removes **0.005 %–0.21 %** of an aldehyde
  during any real thermal step — but over **weeks at ambient in a high-protein matrix it is
  real and sizeable** (a C10 aldehyde loses 14.5–20 % of dose in 28 days at 25 °C, and
  hexanal's half-life there is measured in weeks). If your question is about shelf life
  rather than about the cook, this model answers the wrong question, and it will not warn
  you, because every number it prints is an end-of-cook number.

- **Anything the model refuses.** Do not work around a refusal by substituting a nearby
  compound. Norfuraneol is not DMHF; that substitution is exactly what the refusal exists to
  prevent.

### The single most useful habit

Run the comparison, read the resolved ratios, ignore the absolutes, and book the experiment
that `rank-experiments` puts at the top. The model is a well-instrumented way of being wrong
in a located, quantified way, and that is worth more than a confident number.

---

*Authoritative accuracy figures live in the generated model card in the
[README](../README.md#when-to-trust-the-predictions), which is regenerated from the artifacts.
If this page and the model card ever disagree, the model card is right and this page has
drifted.*
