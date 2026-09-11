# Mittelmaier, Fünfrocken, Fenn, Berlich & Pischetsrieder 2011 — EXTRACTION (3-deoxygalactosone in peritoneal dialysis fluids; and the fed 3-deoxyglucosone experiment that measures the step this model loses)

**Source on disk:** `data/articles/mittelmaier2010.pdf` (0.25 MB; downloaded 2026-09-11 at this
repository's request — it was the top item of the 2026-09-11 reading list). Read 2026-09-11 via
`pdftotext -layout`. Wave B37.

| field | value |
|---|---|
| Title | "3-Deoxygalactosone, a new glucose degradation product in peritoneal dialysis fluids: identification, quantification by HPLC/DAD/MSMS and its pathway of formation" |
| Venue | Analytical and Bioanalytical Chemistry (2011) 399:1689–1697 |
| DOI | 10.1007/s00216-010-4456-3 |
| Group | Friedrich-Alexander-Universität Erlangen-Nürnberg (Pischetsrieder), food chemistry |
| Systems | commercial peritoneal-dialysis fluids; **and a fed-intermediate experiment in a glucose-free PD model at pH 5** |
| What is measured | 3-deoxyglucosone (3-DG), 3-deoxygalactosone (3-DGal) and 3,4-dideoxyglucosone-3-ene (3,4-DGE) as their o-phenylenediamine quinoxalines, by HPLC/DAD with HPLC/DAD/MSMS confirmation, against a five-point calibration with 2,3-dimethylquinoxaline as internal standard |

## 1. Why this paper matters to this model, in one paragraph

The repository's largest single miss is the step **3-deoxyglucosone → 3,4-dideoxyglucosone-3-ene**.
Wave B34 measured it 32× too slow against Leitzen 2021 (autoclaved glucose, 121 °C); the B36
correction measured it 7–10× too slow against Zhang 2021's within-study ratio (aqueous glucose,
90–110 °C). Both of those are pots charged with **glucose**, so the 3-DG → 3,4-DGE step is inferred
from a chain. `docs/guides/EXPERIMENTS.md` therefore asks for "the rate of 3-deoxyglucosone →
3,4-dideoxyglucosone in water, at 100–140 °C and pH 4–7". **This paper runs exactly that experiment:
it charges PURE 3-DG, in water at pH 5, at 120 °C, and follows 3,4-DGE against time.** It is the fed
experiment, not the inferred one.

## 2. The fed experiment, verbatim (sec. "Reaction pathways of 3-DGal formation")

> "To investigate the pathway of 3-DGal formation during sterilization, 3-DG, 3-DGal, or 3,4-DGE
> (∼200 μM each) were dissolved in a PD model without glucose (pH 5). An aliquot of 0.5 mL each of
> the solutions was then heated in gas-tight headspace vials at 120 °C for 0, 10, 20, 30, 60, 90, and
> 120 min. After cooling the samples on ice, quantification was carried out as described above."

So: **~200 µM of ONE pure dicarbonyl, 0.5 mL, sealed gas-tight vial, pH 5, 120 °C, seven time
points to 2 h.** The medium is "a PD model without glucose", i.e. the salts of a conventional
peritoneal-dialysis fluid (the paper's "PD buffer") with the glucose left out. **The paper does not
print the salt composition of that model in this section**, and this dossier does not supply one: a
conventional PD fluid is a lactate-buffered saline, but *likely* is not *stated*, and the
repository's rule is that the medium is recorded as declared or not at all.

## 3. The numbers, verbatim from the Results (the time courses themselves are Fig. 5 — FIGURE ONLY)

| quantity | printed value |
|---|---|
| fed **3-DG** → 3,4-DGE maximum | **26.7 µM at 30 min** |
| fed **3-DGal** → 3,4-DGE maximum | 46.2 µM at 20 min |
| fed **3,4-DGE** → 3-DG maximum | 26.9 µM at 30 min |
| fed **3,4-DGE** → 3-DGal maximum | 37.9 µM at 20 min |
| diastereomer interconversion, fed 3-DG | "the percentage of the respective other diastereomer increased from initially less than 1% to **26% (3-DG)** ... after 60 min of heating" |
| diastereomer interconversion, fed 3-DGal | "... and to **48% (3-DGal)** after 60 min of heating" |

Verbatim on the decay: *"A rapid decline of the overall GDP concentration was observed in all three
experiments. Interestingly, the decrease was much slower when 3,4-DGE was heated compared with 3-DG
and 3-DGal."* And: *"The results confirm that 3-DG, 3,4-DGE, and 3-DGal can be converted into each
other by reversible reactions."*

**Figure 5 carries the full time courses** (panels a, b, c; triangles 3-DG, inverted triangles
3-DGal, circles 3,4-DGE, squares the sum). Only the six maxima above are printed as numbers. The
individual time points are a figure read and are NOT transcribed here.

## 4. What the repository may take, and the three reasons it is not a drop-in rate constant

The transferable quantity is a **fed yield**: from ~200 µM 3-DG at 120 °C and pH 5, the 3,4-DGE
pool peaks at 26.7 µM, i.e. **13.4 % of the charge, at 30 minutes**. Under the owner's rule
(*rates, activation energies, fed yields and within-study ratios FIT; end-of-cook levels VALIDATE*)
a fed yield is fit evidence. Three things must be carried with it or the number will be misused:

1. **The step is REVERSIBLE and the paper proves it.** Feeding 3,4-DGE regenerates 3-DG (26.9 µM at
   30 min). The trunk's `r_tdg_ddg` is one-way. A peak height under a reversible mechanism is not the
   same observable as a peak height under an irreversible one, and fitting the one-way constant to
   this peak without saying so would be transferring a number across a mechanism change.
2. **A peak height constrains a RATIO, not a rate.** The height of an intermediate's maximum is set
   by the formation constant *and* the removal constant together; the time of the maximum carries the
   rest. 26.7 µM at 30 min is therefore a joint constraint on the entry to 3,4-DGE and its exit, and
   the exit in this pot includes a route the model does not have at all (below).
3. **3-deoxygalactosone is not a species in this model.** It is a real sink here — 26 % of the fed
   3-DG has become its C4 epimer by 60 min — and the model has nowhere to put it. Any fit that
   ignores it will push that flux into whichever of the model's channels is nearest, which is exactly
   the failure mode the repository calls a lane conflict.

## 5. The commercial-fluid numbers (context, not a benchmark)

The paper's main result is that 3-DGal is a major, previously unrecognised glucose degradation
product of heat-sterilised PD fluids, formed from 3-DG via 3,4-DGE (Fig. 6 is the mechanism: the
reversible hydration/dehydration of 3-DG ⇌ 3,4-DGE ⇌ 3-DGal). Conventional fluids are sterilised at
**120 °C for approximately 4 h in a buffer at pH 4.5**; low-GDP fluids sterilise the glucose
separately at low pH in a two-compartment bag. Those product surveys are not a cook this repository
can charge and are not taken.

## 6. Verdict

**The most directly useful paper to arrive for this model.** It gives the one experiment the
experiments guide asks for by name, in water, at a cooking temperature, with the intermediate
quantified. It is FIT-class evidence for the `k_tdg_ddg` limb, and it comes with a mechanism finding
(reversibility plus an epimer sink) that is itself a structural question for the trunk. Neither has
been acted on in wave B37: both are pre-registered questions for the wave that follows it.
