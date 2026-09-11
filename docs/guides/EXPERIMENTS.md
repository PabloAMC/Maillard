# What to measure next, and why it would matter

*Written for someone deciding where to spend bench time. Every number in it comes from the model's
own records; nothing here is a guess about what would be interesting.*

## How to read this

A measurement is worth making to this model for one of three reasons, and they are not the same
reason. Knowing which one you are buying tells you what you get back.

**1. It identifies a coordinate the data does not pin.** The model has ten fitted constants that its
own evidence cannot locate: the fit can slide them a long way without the cost changing, so their
values are artefacts of where the search happened to stop. The command `maillard wishlist` lists
them. A measurement here converts a band artefact into a fitted value.

**2. It lifts a refusal.** The model refuses to answer eighteen questions it is asked. Every refusal
names what is missing. Some are missing a number; some are missing a route; one, as of today, is
missing neither and is refused because the lane that would carry it does not drive that row. Read the
refusal before designing the experiment, because two of them cannot be lifted by measuring the thing
they appear to be about.

**3. It is a second laboratory.** Almost every constant in this model comes from one group. Where a
second group has been read, the results have been sobering: three trunk constants agree inside a
factor of two, two declared decisions were refuted outright, and one constant is out by a factor of
466. A replication is worth as much as a new number and costs less.

**One warning about the word "unlocks".** When the model says a measurement would unlock a
prediction, it means the observable sits downstream of that step, so a measured rate would replace a
band artefact with a fitted value. It does **not** mean the answer would then be right. Whether it
lands within threefold of a measurement is what the next pre-registered wave finds out.

---

## The experiments, in the order I would fund them

### 1. Where the thiols go: one pot, two temperatures, fed and unfed in parallel

**The question.** The model removes the meaty thiols far faster than any real pot does, at 100 °C and
at 140 °C alike. One removal step with one temperature dependence cannot fit both. Five different
structures have now been built for that step and all five refused — the fifth without a fit, because
its decisive test turned out to be unreachable by construction (see below).

**Why it cannot be answered by reading.** Four laboratories' papers have been read against this and
none measures a removal rate on a fed thiol with nothing forming. The one relevant storage study is a
declared hold-out, so fitting to it is not allowed.

**The protocol.** Ribose 100 mmol/L and cysteine 33 mmol/L in 0.5 mol/L phosphate at pH 5, in 20 mL
vials with 5 mL of liquid, and **write the volumes down** — headspace volume is an input this model
needs and almost no published pot reports it. Two temperatures: 100 °C from 0.5 to 12 hours, and
140 °C from 5 to 120 minutes. One vial per time point, three replicates. Beside it, the same buffer
with 2-methyl-3-furanthiol alone and 2-furfurylthiol alone at 1 mg/L on the same grid, so removal is
measured with nothing forming. Measure the two thiols by stable-isotope dilution **and their
disulfides in the same run**, plus residual sugar, residual cysteine and the final pH. One arm at
100 °C under nitrogen.

**Add one arm the earlier version of this list did not have.** A fifth set of vials at 140 °C with
the fed thiol plus 2,3-pentanedione at 10 mmol/L. That single addition tests the mechanism the fifth
structure rested on, which is that the pot's own diketones are what oxidise the thiol to its
disulfide. Four papers were fetched to find a rate for that step and not one supplies one: every
constant they print belongs to a competing reaction that makes an adduct instead. So this arm is
currently the only way to get the number.

**What the fifth structure found before it was fitted, and why this arm now matters more.** Building
the diketone oxidant into the model showed the disulfide shortfall is two different problems wearing
one name. In the pots where no oxidant was ever charged, supplying one from the diketones works — but
only if every mercaptoketone-forming event oxidises a thiol, the physical ceiling. In the pots that
already carry air, the model uses under one per cent of its oxidant and still makes ten times too
little disulfide: there the shortfall is the dimerisation *rate*, and that rate cannot be raised
because the same step would then destroy a pure thiol in buffer, which Kumazawa measured surviving.
So the pure-thiol arm and the pure-thiol-plus-diketone arm above are no longer a side check; together
they are the measurement that separates a rate the model has wrong from an oxidant it lacks.

**What the model predicts today.** The fed thiol decays to essentially zero, and the reacting pot
peaks after about an hour and then falls. Both are almost certainly wrong.

**What each outcome decides.** If the fed thiol levels off and its disulfide accounts for the
difference, the removal step is reversible and the disulfide branch is real. If the thiol levels off
and the disulfide does not account for it, the sink is something else and the fourth structure is
still missing. If the diketone arm makes markedly more disulfide than the arm without it, the
oxidant hypothesis is right and the constant falls straight out of the comparison.

**Cost.** About 130 vials and two weeks of gas chromatography.

### 2. A fed-intermediate ladder: six unidentified constants in one design

**The question.** Ten fitted constants are unidentified. Six of them ask for exactly the same kind of
measurement, so they can be bought together rather than one at a time.

**The protocol.** Heat each of the following alone in the lane's own buffer at its reference
temperature of 145 °C, and quantify the named product against time by stable-isotope dilution where
a labelled standard exists: the Amadori compound alone, measuring the deoxypentosone; furfural alone,
measuring its loss; furfural with hydrogen sulfide, measuring 2-furfurylthiol; each thiol's disulfide
alone, measuring its loss; glucose alone, measuring hydroxyacetaldehyde and methylglyoxal. Six
sub-experiments, one buffer, one temperature, one analytical method.

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
corpus — six of them — computes odour activity in a plant matrix using a **water** threshold, and
most take it from the same database. Two papers were fetched specifically to lift this refusal. One
turned out to compare soymilk against soymilk with no water leg at all. The other is a genuine paired
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
explain — a cap computed from one compound in beef and one dairy protein, and exceeded for the first
time this week on a real hold-out row — transfers to plant protein at all.

### 4. A flavour binding constant on a plant protein, in water, hot

**The question.** The model's binding constants for plant protein are measured at 37 °C, a mouth
temperature. It says plainly that nothing licenses them at 90 or 140 °C.

**What the reading did and did not settle.** The only flavour-affinity measurement in the corpus above
60 °C was fetched this week. It runs at 80, 90 and 100 °C, and it is **dry**: gas-solid chromatography
on a powder, no water anywhere, and its own authors write that the binding order will change in
solution. It checks the model's chain-length slope, agreeing to within 1.26-fold, and it supplies no
constant this layer can use. So the gap is real and reading has now exhausted itself against it.

**The protocol.** Phase-ratio variation on a pea protein isolate at a stated loading, against hexanal
and one alkenal, at 40, 70 and 90 °C, with a water leg measured in the same run on the same
instrument. The water leg is the part that matters: the model's stored form is a ratio precisely
because the absolute headspace scale is untrustworthy, and a paper that prints only a matrix leg
cannot be used no matter how carefully it was done. One paper fetched this week failed exactly there.

**What it decides.** Whether binding weakens with heat, as the dry measurement suggests, or
strengthens, as two aqueous studies on other proteins report. The model currently charges its binding
sites once at the start of a cook and does not change them, and has no evidence either way.

### 5. The melanoidin's carbon-to-nitrogen ratio, as a series

**The question.** The brown polymer's composition is the model's one direct check on whether its
browning arithmetic is right. The model assumes a fixed repeat unit of eight carbons per nitrogen.

**What the reading settled, and it is not what was expected.** Five laboratories now bracket this and
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

## Cheaper than a laboratory

Three things would help and cost an email, and the first is the single cheapest way to make this
model answer more questions than any experiment on the list above.

**Run your blank.** If you are measuring volatiles in a heated plant protein, measure the SAME
material unheated, in the same run, and print that column beside the heated one. It costs one extra
vial. What it is worth, measured rather than argued:

A pea beverage paper printed both columns. Between 36 % and 42 % of every volatile level it reports
for the heated product was already in the beverage before any heat — hexanal 331 of 782 µg/L,
2-pentylfuran 59.4 of 163, nonanal 8.24 of 24.0. Scored against the totals, this model missed those
three by 34x, 32x and 3.3x. Told what the beverage started with, it lands at **2.2x, 2.5x and 1.5x**,
all three inside the threefold band, on rows no fit has read. Nothing in the model changed. It simply
stopped being charged for raw material it never had to make.

The other side of the same coin is four pots this model now REFUSES. They are unheated protein
powders and flours whose measured hexanal is entirely what the ingredient arrived with, and no source
on file prints anything prior to compare against. Asked to form those levels from zero in ten
minutes, the model missed by factors of 3 357, 6 078, 3 717 and 33 392. Those are not chemistry
failures, and reporting them as such would be dishonest, so the rows are refused — and the refusal
prints the cure, which is this blank. **Supply an unheated column for those materials and seven
refused rows become answerable.** Nothing else on this page unlocks that many for that little.

The general rule is not specific to this model: **a formation measurement without its own blank
cannot be told apart from a storage measurement.** If your material sat in a warehouse for six
months, some of what you are about to attribute to your process was in the bag when you opened it.

**One aqueous rate, and it is now a different one than this page said yesterday.** The entry from a
sugar to 3-deoxyglucosone was the suspect. It is not: measured against a paper that turned out to be
on this disk already, the model's 3-deoxyglucosone in an amine-free pot is right to **11 %** and its
methylglyoxal to **28 %**. What is wrong is the step immediately after — 3-deoxyglucosone to
3,4-dideoxyglucosone, **32× too slow**, with hydroxymethylfurfural 12× low downstream of it.

So the ask is narrow and specific: **the rate of 3-deoxyglucosone → 3,4-dideoxyglucosone in water,
at 100–140 °C and pH 4–7.** The model carries a value for it, and that value's interval was
deliberately *not* widened because two laboratories agree on it to within 1.5×. Both of them measured
it dry — a freeze-dried glucose melt and a roasted hazelnut. Two laboratories agreeing with each
other in the same matrix is not a check; it is the same experiment twice. One aqueous determination
would settle it, and it is the cheapest number on this page.

A candidate exists and is open access, not yet on this disk: Zhang, Sun, Pu, Zhang, Sun & Zhao 2021,
*Food Science & Nutrition* 9:290–302 (doi 10.1002/fsn3.1995) heat 0.3 M glucose **in water, alone and
with glutamic acid, at 90–110 °C for 0–6 h** and quantify 3-deoxyglucosone, 1-deoxyglucosone,
**3,4-dideoxyglucosone**, glucosone, glyoxal, methylglyoxal and diacetyl, with zero-order formation
constants printed per temperature. It measures the very intermediate this model loses; whether its
time courses support a 3-DG → 3,4-DGE rate is the first thing to check once it is read.

**Numbers behind figures.** Two published data sets that bear directly on the thiol problem exist only
as figures. Their authors' underlying numbers would give that step its first data from a third
laboratory.

**One supplementary table.** The melanoidin series above exists in part already: a paper measures ten
temperatures and prints only the two endpoints, with the rest in a supplementary file that is not
publicly posted. Asking for it would supply most of experiment 5 for nothing.

---

## What not to spend money on, and why

**Do not measure a per-amino-acid Strecker rate constant in water.** Four laboratories have now been
read against this gap and every one printed a barrier or a yield instead of a rate. The likeliest
reading is that nobody has measured it because the measurement is hard, not because it was missed.
That makes it a research project, not a gap to fill.

**Do not measure 2-pentylfuran's branch fraction.** It was measured in 1981 and entered the model.
This entry previously said the compound was still refused because the rows asking for it are scored
against a hexanal the lipid lane does not produce. **That was wrong and is retracted.** The lane's
answers looked a hundred thousand times too small because the engine had no unit-conversion entry for
the new species and was reporting it in mmol/L instead of µg/L. With that fixed the compound is
answered, and in the one pot that declares its own starting state it lands at 2.5x. The branch
fraction is not the gap. Where the alkylfuran still misses badly — 366x in an extruded soy row — the
gap is the size of that matrix's hydroperoxide pool, which is a declared input nobody has measured
for it, and a peroxide value would close it.

**Do not add more directional claims to the moisture axis.** That axis reads "do not use", and more
agreeing claims will not lift it, because only one of the model's four lanes carries a water-activity
term at all. The others refuse the comparison structurally. The axis needs a moisture-dependent step
measured and fitted on those lanes, not more observations.

---

## Feeding results back

The model scores your measurements without being changed by them, and calibrates to your laboratory
as a separate file you apply. One laboratory's four-temperature series, two pots fitted and two held
out, moved the held-out pots from a hundredfold error to within threefold, with the shipped model
untouched. See the [quick start](QUICKSTART.md) for `score` and `calibrate`, and
[the wishlist](../../results/validation/data_wishlist.md) for the machine-generated list this guide
is written from.
