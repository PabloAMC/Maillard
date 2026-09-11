# Luo, Tian, Li, Zhang, Bi, Fu & Jing 2024 — EXTRACTION (a mechanism review of volatile sulfur compound formation; no kinetics, one corroboration that matters)

**Source on disk:** `data/articles/Luo2024.pdf` (3.0 MB; downloaded 2026-09-11). Read 2026-09-11 via
`pdftotext -layout`. Wave B43.

| field | value |
|---|---|
| Title | "Mechanisms underlying the formation of main volatile odor sulfur compounds in foods during thermal processing" |
| Venue | Comprehensive Reviews in Food Science and Food Safety 23(4), 2024 |
| DOI | 10.1111/1541-4337.13389 — read from the printed header |
| Group | Henan Agricultural University and Beijing Technology and Business University |
| Type | **COMPREHENSIVE REVIEW.** Eighteen mechanism figures, no experiments of its own |

## 1. The verdict first

**This paper contains no rate constant, no activation energy, no time course and no temperature
series.** Searched for all four; nothing. It cannot move a number in this model and nothing here
is fitted. What it is good for is the layer this repository built for exactly this kind of source:
`data/lit/reaction_rules.yml`, the cited reaction rules run over the model's species to list steps
the literature draws that the engine does not have.

## 2. The corroboration that matters, verbatim

> "The anaerobic oxidation of thiols mainly involves **dicarbonyl compounds as oxidants**. In this
> process, thiols are converted into disulfides, whereas carbonyl groups are reduced to hydroxyl
> groups."

That sentence is the general form of the structure **wave B27 built and then gated**
(`sulfur.ch_redox_mp3p`, {NF, H₂S} → {MP3P, OX} on `k_redox_mp3p`, from Whitfield 1999's Figure 6).
B27 refused it on evidence — in the pots that carry air the dimer step uses under 1 % of its
oxidant, and the fed pots need the coupling at its physical ceiling — and left it inert. **A review
independently draws the same mechanism.** That changes the step's status from one laboratory's
figure to a reviewed mechanism and changes nothing about its rate, so B27's refusal stands
unaltered: the wave was gated on the size of the flux, not on whether the arrow exists.

## 3. The other thiol pathways it draws (Figure 7), and where each already sits

| pathway, as the review draws it | this repository |
|---|---|
| aerobic oxidation, thiol + O₂ → disulfide, "occurs rapidly at room temperature" | rule `R12_thiol_oxidation_to_disulfide`; the lane's oxygen pools are B11's `OX`/`OXR`/`OXV` |
| **anaerobic oxidation, dicarbonyls as the oxidant** | B27's `ch_redox_mp3p`, built inert and gated (above) |
| transition-metal catalysis: Cu²⁺ and Fe³⁺/Fe²⁺ accelerate the oxidation; a Fe²⁺ + H₂O₂ + O₂ Fenton system "rapidly oxidizes thiols into disulfide and even sulfoxide" | already on disk as **kinetics** and already refused: `bagiyan2004_extraction.md` gives Cu-catalysed initial rates at one temperature (20 °C) whose own authors disclaim a quantitative rate law |
| thiol + aldehyde by Michael addition → methylthio-alcohol, which may be substituted again to a dialkyl disulfide | rule `R13_thiol_michael_addition` (the second substitution is not modelled) |
| thiol + alcohol or alkene → alkyl sulfide or disulfide by nucleophilic substitution | not modelled; no rate in the review |

## 4. What it changes here

One line in `reaction_rules.yml`: the disulfide rule's condition string now names the two oxidants
this review distinguishes, with Luo as the source for the anaerobic one, and B27's declaration cites
it as corroboration. **No constant moves, no rule is added, no refusal is lifted.** A review is
evidence that an arrow is drawn, never that it is fast.

## 5. And one thing it makes sharper for the experiment on the list

The metal row above is not idle. This repository's own buffer note on the Yiltirak pots records that
their buffer was made in **tap water**, "so trace-metal catalysis is uncontrolled" — and the sulfur
lane's cross-laboratory scatter is the model's largest single problem. Luo's review says trace copper
and iron change the thiol oxidation rate substantially. That is a candidate explanation for scatter
that no barrier can absorb, and it is cheap to test: **the thiol experiment in `EXPERIMENTS.md` should
carry a chelator arm.** Recorded there.
