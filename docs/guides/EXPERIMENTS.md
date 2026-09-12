# What to measure next, and why it would matter

*Written for someone deciding where to spend bench time. Every number in it comes from the model's
own records; nothing here is a guess about what would be interesting. The history of how each entry
came to say what it says lives in the pre-registrations under `results/validation/`, not here.*

## How to read this

A measurement is worth making to this model for one of three reasons, and they are not the same
reason. Knowing which one you are buying tells you what you get back.

**1. It identifies a coordinate the data does not pin.** The model has ten fitted constants that its
own evidence cannot locate: the fit can slide them a long way without the cost changing, so their
values are artefacts of where the search happened to stop. The command `maillard wishlist` lists
them. A measurement here converts a band artefact into a fitted value.

**2. It lifts a refusal.** The model refuses to answer eighteen questions it is asked. Every refusal
names what is missing. Some are missing a number; some are missing a route; one is missing neither
and is refused because the lane that would carry it does not drive that row. Read the refusal before
designing the experiment, because two of them cannot be lifted by measuring the thing they appear to
be about.

**3. It is a second laboratory.** Almost every constant in this model comes from one group. Where a
second group has been read, the results have been sobering: three trunk constants agree inside a
factor of two, two declared decisions were refuted outright, and one constant is out by a factor of
466. A replication is worth as much as a new number and costs less.

**One warning about the word "unlocks".** When the model says a measurement would unlock a
prediction, it means the observable sits downstream of that step, so a measured rate would replace a
band artefact with a fitted value. It does **not** mean the answer would then be right. Whether it
lands within threefold of a measurement is what the next pre-registered wave finds out.

**One rule that applies to every experiment below.** The only fit on this branch that pinned every
constant it touched was built on a paper that fed a pure intermediate and followed it against time.
Fits scored on end-of-cook levels through the whole network leave most of their constants unseen by
any row, and the identifiability audit says the same from the Fisher information: it is the design
that pins, not the number of rows. So every protocol here feeds something pure and follows it. If you
adapt one, keep that shape.

---

## The experiments, in the order I would fund them

### 1. Where the thiols go: one pot, two temperatures, fed and unfed in parallel

**The question.** The model is wrong about the meaty thiols in two directions at once, and both have
to be measured in the same design or neither can be told from the other.

A *fed* thiol, put in buffer with nothing forming, decays in the model to essentially zero, far
faster than any real pot loses it. Five different structures have been built for that removal step
and all five were refused. The identifiability audit then found the step cannot be moved by any refit:
its barrier and both dimerisation rates sit on their ceilings together, and raising the rate would
destroy a pure thiol in buffer that Kumazawa measured surviving.

A *reacting* pot of ribose and cysteine, on the other hand, comes out too high. Across nineteen scored
rows from eleven separate pots the model **reads high on both meaty thiols, not low**, roughly
thirtyfold on the roasted one and sixfold on the meaty one, on the same side in four rows out of
five. A sink that is too weak would produce that, but the sink is pinned and cannot be made stronger
on the present evidence. **The formation side is now equally implicated, and no wave has examined
it.** The over-prediction is largest at the lowest temperatures and shrinks as the pot gets hotter.

The two findings fit together. A probe of cysteine alone in buffer showed the model keeping 99.5 % of
it after five minutes at 95 °C, and 83 % after three hours, where a real pot carrying trace copper has
lost nearly all of it within five minutes of a gentle ramp to 60 °C. Cysteine the model fails to lose
is cysteine still available to make thiols. So one missing removal step upstream shows up downstream
as over-production, and it matters most where thermal chemistry is slowest, which is exactly where
the offset is largest. What is missing is a catalytic channel whose rate depends on a catalyst the
model does not carry. No adjustment of a thermal constant can imitate one.

**Why it cannot be answered by reading.** Six laboratories' papers have been read against this. None
measures a removal rate on a fed thiol with nothing forming in a defined buffer; the closest starts
every run from a thioacetate and its loss rate scales with the dose of a crude enzyme. The one
relevant storage study is a declared hold-out, so fitting to it is not allowed. And the one paper that
measured trace-metal catalysis of a thiol did so on a thioether and never assayed the thiol in its
chelator arm.

**The protocol.** One buffer, two temperatures, five arms, one analytical method.

- **Buffer and pots.** Ribose 100 mmol/L and cysteine 33 mmol/L in 0.5 mol/L phosphate at pH 5, in
  20 mL vials with 5 mL of liquid. **Write the volumes down.** Headspace volume is an input this model
  needs and almost no published pot reports it. State the water's provenance; a buffer in this
  model's own corpus was made in tap water and its trace-metal content is unknown.
- **Time grid.** 100 °C from 0.5 to 12 hours, and 140 °C from 5 to 120 minutes. One vial per time
  point, three replicates.
- **Arm A, the reacting pot**, as above.
- **Arm B, the fed thiols.** The same buffer with 2-methyl-3-furanthiol alone and 2-furfurylthiol alone
  at 1 mg/L, on the same grid, so removal is measured with nothing forming.
- **Arm C, oxygen.** Arm B repeated at 100 °C under nitrogen.
- **Arm D, the chelator.** Arms A and B repeated with a metal chelator, and **the chelator must be in
  molar excess over the thiol, not over the metal.** A chelator at one-to-one with the copper and
  short of the thiol makes the oxidation worse, not better; the loss it was meant to suppress roughly
  doubles, because chelation only delays the metal's transfer to the thiol until a higher temperature.
  Only a genuine excess abolished the effect. Do not assume the arm without chelator is clean because
  the water was: in ultrapure water every transition metal sat below detection until cysteine was
  added, and the cysteine itself, at the highest purity grade sold, carried enough copper to raise the
  solution to a quarter of a micromolar and drive a resolved oxidation at 95 °C. The reagent is a metal
  source.
- **Arm E, the diketone.** Arm B at 140 °C with 2,3-pentanedione at 10 mmol/L added. This tests
  whether the pot's own diketones are what oxidise the thiol to its disulfide. Four papers were fetched
  for that rate and every constant they print belongs to a competing reaction that makes an adduct
  instead, so this arm is currently the only way to get the number.

**What to measure, in order of importance.**

1. **Residual cysteine against time, on arm D and its counterpart without chelator.** This is the
   pivot, not a housekeeping check. The difference between those two cysteine curves is the size of
   the missing channel.
2. **Both thiols and both disulfides in the same run**, by stable-isotope dilution. Without the
   disulfide in the same run the mass balance cannot be closed, and a thiol that disappears cannot be
   attributed to dimerisation rather than to volatilisation or something else.
3. Residual sugar and the final pH.

**What each outcome decides.**

- If cysteine disappears faster without the chelator and the thiols come down with it, the
  over-production is explained by a missing catalytic removal of the precursor, and the sink was never
  the main story.
- If the fed thiol levels off and its disulfide accounts for the difference, the removal step is
  reversible and the disulfide branch is real.
- If the fed thiol levels off and the disulfide does not account for it, the sink is something else
  and its structure is still missing.
- If arm E makes markedly more disulfide than arm B, the diketone oxidant is real and its constant
  falls straight out of the comparison. Building that oxidant into the model showed the disulfide
  shortfall is two problems wearing one name: in pots that were never charged with oxygen the
  diketone supply works, but only at its physical ceiling; in pots that carry air the model uses under
  one per cent of its oxidant and still makes ten times too little disulfide, so there the shortfall is
  the rate, and the rate is pinned. Arms B and E together are what separate a rate the model has wrong
  from an oxidant it lacks.

**Cost.** About 130 vials and two weeks of gas chromatography.

### 2. A fed-intermediate ladder: seven constants in one design

**The question.** Ten fitted constants are unidentified. Six of them ask for exactly the same kind of
measurement, so they can be bought together rather than one at a time. A seventh measurement of the
same shape, on a different lane, would give three recently added steps their own temperature
dependence.

**The protocol.** Heat each of the following alone in the lane's own buffer and quantify the named
product against time by stable-isotope dilution where a labelled standard exists.

At the pentose lane's reference temperature of 145 °C:

- the Amadori compound alone, measuring the deoxypentosone;
- furfural alone, measuring its loss;
- furfural with hydrogen sulfide, measuring 2-furfurylthiol;
- each thiol's disulfide alone, measuring its loss;
- glucose alone, measuring hydroxyacetaldehyde and methylglyoxal.

And on the sugar lane, at 100 and 140 °C:

- pure 3-deoxyglucosone in water at pH 5, following 3,4-dideoxyglucosone and 3-deoxygalactosone.
  This experiment has been run once, at 120 °C, and it changed the model: the dehydration runs both
  ways, the enone hydrates to either sugar epimer, and a quarter of the fed compound is the galactose
  epimer after an hour. Three steps were added and fitted from it, all pinned to a tenth of a decade,
  and they share one declared barrier because their source has one temperature. Two more temperatures
  on the same pots would give each its own, and cost one afternoon. What is left on that limb after
  the fit is a residual 6.6× on 3,4-dideoxyglucosone at 121 °C and pH 4.4, which is an extrapolation
  of a pH term measured between 5.5 and 6.8.

**What it unlocks, in the model's own words.** The glucose one alone would return eight panel rows
that are currently answered but declared not evaluable, and would give absolute answers from any
glucose-containing charge instead of a band artefact. The others each remove one sliding coordinate.

**Why it is second and not first.** It makes the model more precise where it is already roughly
right. Experiment 1 decides whether a whole structural claim is right or wrong. Precision after
correctness.

### 3. Paired odour thresholds in a plant protein

**The question.** The model refuses every matrix-corrected odour threshold. It can tell you a
compound's concentration and its threshold in water, and it refuses to tell you what that threshold
becomes in a pea or soy dispersion.

**Why it cannot be answered by reading, and this was checked hard.** Every plant-protein paper in the
corpus, six of them, computes odour activity in a plant matrix using a **water** threshold, and most
take it from the same database. Two papers were fetched specifically to lift this refusal. One turned
out to compare soymilk against soymilk with no water leg at all. The other is a genuine paired
threshold and it is in whey, not a plant protein, where the shift is somewhere between eightyfold and
eleven thousandfold depending on which water baseline you accept. There is no plant-matrix paired
threshold anywhere in reach. **This refusal needs a sensory panel and nothing else will do.**

**The protocol.** Six to eight compounds spanning the classes the model emits: an n-alkanal
(hexanal), a branched aldehyde (3-methylbutanal), a methyl ketone, an alkylpyrazine, a thiol
(2-furfurylthiol), a disulfide (dimethyl trisulfide), and a methoxypyrazine. Detection thresholds by
the ascending three-alternative forced-choice method, orthonasal, with the standard correction for
guessing, in water and in a 3 % w/w pea protein dispersion at the same pH, on the same panel in the
same sessions. Report the panel size and the criterion, because the corpus is full of thresholds that
cannot be compared for want of them.

**What it decides.** Whether the model's cap on how much of a matrix shift reversible binding can
explain, a cap computed from one compound in beef and one dairy protein and already exceeded on a real
hold-out row, transfers to plant protein at all.

### 4. A flavour binding constant on a plant protein, in water, hot

**The question.** Two things, and the second is larger than the first.

The model's binding constants for plant protein are measured at 37 °C, a mouth temperature. It says
plainly that nothing licenses them at 90 or 140 °C.

And the model holds exactly one way for an aldehyde to attach to a protein: a permanent chemical bond
to the lysines, with a rate measured on dairy proteins at room temperature. That is not the channel a
headspace measurement actually sees. Acidifying a soy isolate and holding it at 95 °C for five minutes
nearly **tripled** the hexanal in the headspace, and in the storage globulin alone it rose more than
fivefold: heat *releasing* aldehyde, not consuming it. Run against that same pot, the model moves
hexanal by about one part in ten thousand, downward. The sign is wrong and the size is out by roughly
four orders of magnitude. The permanent bond is real and its measured rate is not in dispute. Sitting
on top of it is a much larger, **reversible** association, weak, physical, undone by acid and heat
together, that governs how much aldehyde ever reaches the headspace. The model does not have it.

**Why it cannot be answered by reading.** The only flavour-affinity measurement in the corpus above
60 °C runs at 80, 90 and 100 °C and is **dry**: gas-solid chromatography on a powder, no water
anywhere, and its own authors write that the binding order will change in solution. It checks the
model's chain-length slope, agreeing to within 1.26-fold, and supplies no constant this layer can
use. The soy study that showed the release measures headspace only, never the total, so it cannot
give a constant either. Reading has exhausted itself against this gap twice.

**The protocol.** Phase-ratio variation on a pea protein isolate at a stated loading, against hexanal
and one alkenal, at 40, 70 and 90 °C, with three requirements.

- **A water leg measured in the same run on the same instrument.** The model's stored form is a
  ratio precisely because the absolute headspace scale is untrustworthy, and a paper that prints only
  a matrix leg cannot be used no matter how carefully it was done. One paper fetched for this gap
  failed exactly there.
- **Total hexanal on the same aliquot as the headspace value**, by exhaustive extraction or by
  purging to completion, at both pH 4.5 and pH 7. Without that pairing the experiment returns a
  constant for the wrong channel.
- **Hexanal against time at the three temperatures, in the same suspension.** The same pots answer a
  second question for almost nothing: the fat path's temperature dependence in a hot, wet protein. It
  has been measured twice and the two answers disagree. In bulk seed oil the barrier for making
  hexanal from an existing peroxide pool is about 114 to 122 kilojoules per mole between 130 and
  160 °C; in a moist nut paste between 60 and 130 °C it is 114 when the matrix is dry and falls to
  about 62 once it holds water. The model bridges its room-temperature anchor to cooking temperature
  with a single factor that sits above every wet-matrix measurement at every temperature. Which of the
  two applies to a plant protein in water is unknown, and three papers fetched against exactly that
  question all failed, one on design and two on how they reported.

**What it decides.** Whether binding weakens with heat, as the dry measurement suggests, or
strengthens, as two aqueous studies on other proteins report; how large the reversible channel is and
whether acid opens it; and which of the two published barriers governs hexanal in a wet protein.

**What this costs the model until it is done.** In a protein matrix, a headspace measurement reports
the share that escapes and this model reports the whole amount, and the gap between them moves with
the sample's acidity and its heating history. Every hexanal row measured that way is therefore
compared a little unfairly. It would have been easy to widen those rows' tolerances on this reasoning
and let the scores improve; that was refused, because it would raise the headline without
establishing a single new fact about the chemistry. The mismatch is declared instead, and this
measurement is what removes it.

### 5. The melanoidin's carbon-to-nitrogen ratio, as a series

**The question.** The brown polymer's composition is the model's one direct check on whether its
browning arithmetic is right. The model assumes a fixed repeat unit of eight carbons per nitrogen.

**What the reading settled, and it is not what was expected.** Five laboratories bracket this and
they falsify the fixed unit **from both directions at once**: below it at 70 °C, where two thirds of
the amine arrives with its carboxyl removed and contributes one carbon instead of two, and far above
it from 100 °C upward, where the polymer takes up more sugar per amine than one unit allows. A single
fixed unit cannot be wrong in both directions and still be the right structure. What the evidence
supports is a sugar-to-amine ratio that rises with temperature, which the model does not have.

**The protocol.** Glucose and glycine, equimolar, in one buffer at one pH, taken to the same browning
endpoint at 70, 90, 110 and 130 °C. Dialyse identically at every temperature. Elemental analysis, plus
a labelled-carbon or carbon-dioxide measurement so the decarboxylated fraction is measured rather than
inferred. The corpus has this at one temperature in one system and everything else is cross-study.

**What it unlocks.** A temperature-dependent ratio to replace a fixed unit, which is a wave the model
cannot run today for want of exactly this series.

---

## Without a laboratory

Four things cost an email or an extra vial, and the first makes this model answer more questions than
any experiment above.

**Run your blank.** If you are measuring volatiles in a heated plant protein, measure the same
material unheated, in the same run, and print that column beside the heated one. It costs one extra
vial. What it is worth, measured rather than argued:

A pea beverage paper printed both columns. Between 36 % and 42 % of every volatile level it reports
for the heated product was already in the beverage before any heat: hexanal 331 of 782 µg/L,
2-pentylfuran 59.4 of 163, nonanal 8.24 of 24.0. Scored against the totals, this model missed those
three by 34×, 32× and 3.3×. Told what the beverage started with, it lands at **2.2×, 2.5× and 1.5×**,
all three inside the threefold band, on rows no fit has read. Nothing in the model changed. It simply
stopped being charged for raw material it never had to make.

Four pots this model refuses are the other side of the same coin. They are unheated protein powders
and flours whose measured hexanal is entirely what the ingredient arrived with, and no source on file
prints anything prior to compare against. Asked to form those levels from zero in ten minutes, the
model missed by factors of 3 357, 6 078, 3 717 and 33 392. Those are not chemistry failures, so the
rows are refused, and the refusal prints the cure. **Supply an unheated column for those materials and
seven refused rows become answerable.**

A blank corrects the target; it does not cure a process declared as one hold. The extrusion paper
behind this model's only acrylamide process row prints its unextruded control as a bar, about 38 µg/kg
against about 150 in the extrudate, so a quarter of that row's target was in the bag. It is not
declared, because it is read off a figure, and because the row misses by 4 247×, which no blank can
fix: the process is scored on its die zone alone.

The general rule is not specific to this model. **A formation measurement without its own blank
cannot be told apart from a storage measurement.** If your material sat in a warehouse for six months,
some of what you are about to attribute to your process was in the bag when you opened it.

**Ask one author for eleven numbers.** A paper on reconstituted whey protein ran the exact experiment
the fat path needs: hexanal in real concentration units at seven temperatures from 50 to 90 °C in a
wet protein solution at neutral pH, with two holding times. Every value is a bar in a figure. The
article's supplementary tables are an identification list and a set of calibration curves, not the
concentrations, so a download will not do it; the numbers most likely exist only in the authors'
records. Ask the corresponding author for the means and standard deviations behind that figure, and
for an unheated control. Fetch the supplementary file too, because the calibration curves and the
heating-and-cooling profile are both needed for a fit. Even complete, the result would carry a wide
band: there is no unheated baseline, and the protein contains no measurable fat, so the substrate is
unquantified. What it buys is a third laboratory on the one axis where this model currently
interpolates between two that disagree, and what can already be taken from it without reading a bar
height is the shape: across a 40 K span, eight of its eleven treatments are statistically
indistinguishable, the signature of a weak temperature dependence.

**Ask for two supplementary tables that exist.** One holds the concentrations of the two meaty thiols
and all three of their dimers in one pot across five reaction times; the published article prints
only odour-activity ratios. The other holds a plant protein's volatiles at four temperatures with an
unheated blank; the published article prints only the class totals. Both experiments have been done.
Either would answer a question on this page outright.

**Ask for the melanoidin series.** A paper measures experiment 5's ratio at ten temperatures and
prints only the two endpoints, with the rest in a supplementary file that is not publicly posted.
Asking for it would supply most of experiment 5 for nothing.

---

## What not to spend money on, and why

**Do not measure a per-amino-acid Strecker rate constant in water.** Four laboratories have been read
against this gap and every one printed a barrier or a yield instead of a rate. The likeliest reading is
that nobody has measured it because the measurement is hard, not because it was missed. That makes it
a research project, not a gap to fill.

**Do not measure 2-pentylfuran's branch fraction.** It was measured in 1981 and entered the model. An
earlier version of this entry said the compound was still refused for a routing reason; that was
wrong and is retracted. The lane's answers looked a hundred thousand times too small because the
engine had no unit-conversion entry for the species and was reporting it in mmol/L instead of µg/L.
With that fixed the compound is answered, and in the one pot that declares its own starting state it
lands at 2.5×. Where it still misses badly, 366× in an extruded soy row, the gap is the size of that
matrix's hydroperoxide pool, which is a declared input nobody has measured for it. **A peroxide value
on that matrix would close it**, and a peroxide value is the single largest assumption on the whole
fat path.

**Do not add more directional claims to the moisture axis.** That axis reads "do not use", and more
agreeing claims will not lift it, because only one of the model's four lanes carries a water-activity
term at all. The others refuse the comparison structurally. The axis needs a moisture-dependent step
measured and fitted on those lanes, not more observations.

**Do not run a hexanal loss series on an isolate expecting a rate law.** Four papers show hexanal
falling when an isolate is heated, and the reason is now clearer: it is mostly not chemistry. An
enzyme is switched off, a volatile is stripped, an aldehyde is bound. The measurement that separates
the three is the blank above, run with the headspace vented and unvented, not a rate fit.

---

## Feeding results back

The model scores your measurements without being changed by them, and calibrates to your laboratory
as a separate file you apply. One laboratory's four-temperature series, two pots fitted and two held
out, moved the held-out pots from a hundredfold error to within threefold, with the shipped model
untouched. See the [quick start](QUICKSTART.md) for `score` and `calibrate`, and
[the wishlist](../../results/validation/data_wishlist.md) for the machine-generated list this guide
is written from.
